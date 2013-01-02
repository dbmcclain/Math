;; interp.lisp -- Interpolation and table lookup
;;
;; DM/MCFA  04/02
;; -------------------------------------------------

;; A <SLINE-INTERPOLATOR> contains the X, Y, and second derivative table Y2
;; in ascending X order. This object is utilized by SPLINT to provide
;; cubic-spline interpolation in the table of Y(X).
;;
;; It is also used to cache prior lookup values in case the next lookup
;; falls between the same two table nodes. Random lookup is provided in
;; SPLINT by means of binary search for the enclosing nodes.
;;
(defclass <spline-interpolator> ()
  ((xtable   :accessor spl-xtable
             :initarg  :xtable
             :type     (array float (*)))
   (ytable   :accessor spl-ytable
             :initarg  :ytable
             :type     (array float (*)))
   (y2table  :accessor spl-y2table
             :initarg  :y2table
             :type     (array float (*)))
   (nlast    :accessor spl-nlast
             :initarg  :nlast
             :type     fixnum)

   ;; cached entries from previous lookup
   (x0-last  :accessor spl-x0-last
             :initform nil
             :type     (or null float))
   (x1-last  :accessor spl-x1-last
             :initform nil
             :type     (or null float))
   (y0-last  :accessor spl-y0-last
             :initform 0.0
             :type     float)
   (y1-last  :accessor spl-y1-last
             :initform 0.0
             :type     float)
   (y20-last :accessor spl-y20-last
             :initform 0.0
             :type     float)
   (y21-last :accessor spl-y21-last
             :initform 0.0
             :type     float)
   (h-last   :accessor spl-h-last
             :initform 0.0
             :type     float)))


(defun spline (x y &key (kind :natural) yder1 yder2)
  ;; -------------------------------------------------
  ;; This routine should be called only once ahead of
  ;; all the interpolation calls. It creates a triple of
  ;; interpolation vectors suitable for use in SPLINT.
  ;;
  ;; Tables X and Y should be the same sized vectors of numbers.
  ;;
  ;; The X table must be monotonically increasing or decreasing
  ;; in every interval.
  ;;
  ;; We convert them into vectors of floats such that the X table
  ;; is in ascending order.
  ;;
  ;; This routine returns a <SPLINE-INTERPOLATOR> object suitable
  ;; for repeated lookups using routine SPLINT.
  ;;
  ;; Adapted from "Numerical Recipes in C - Second Edition",
  ;; by Press, Teukolsky, Vetterling, and Flannery,
  ;; Cambridge University Press, 1996, Section 3.3 and 3.4, pp. 113-117.
  ;;
  ;; DM/MCFA  04/02
  ;; -------------------------------------------------
  ;; KIND is a keyword indicating the type of cubic spline to provide:
  ;;   :NATURAL sets the endpoint second derivatives to zero
  ;;   :EXTENSIONAL sets the endpoint second derivatives such that the
  ;;                first derivatives at the endpoints agree with the
  ;;                linear slope derived from the pair of nodes at each
  ;;                end.
  ;;   :CLAMPED allows one to specify one or both endpoint slopes,
  ;;            YDER1 and YDER2. If not specified then an extensional slope
  ;;            is adopted. (see description with :EXTENSIONAL).
  ;; -------------------------------------------------
  
  ;; convert to uniform key representation
  (ecase kind
    (:natural)
    (:clamped)
    (:extensional
     (setf kind :clamped
           yder1 nil
           yder2 nil)))

  ;; ensure the claming slopes are float values if specified
  (let ((yder1 (if yder1 (sfloat yder1)))
        (yder2 (if yder2 (sfloat yder2))))

    (declare (type (or null float) yder1 yder2))
  
    ;; ensure we have vectors of floats
    (let* ((x (map 'vector #'single-float x))
           (y (map 'vector #'single-float y)))
      (declare (type (array single-float (*)) x y))
      
      ;; make working vectors
      (let* ((nel   (length x))
             (nlast (1- nel))
             (y2    (make-array nel
                                :element-type 'single-float
                                :initial-element 0.0))
             (u     (make-array nel
                                :element-type 'single-float
                                :initial-element 0.0)))
        
        (declare (type fixnum nel nlast)
                 (type (array single-float (*)) y2 u))
        
        ;; ensure that x and y tables are the same size
        (assert (= nel (length y)))
        
        ;; if handed a reversed vector in x
        ;; then form the forward order x and y vectors
        (destructuring-bind (x y)
            (if (> (aref x nlast)
                   (aref x 0))
                (list x y)
              (list (reverse x) (reverse y)))
          
          (declare (type (array single-float (*)) x y))
          
          ;; ensure that the table is monotonically increasing in x
          (assert (every #'> (subseq x 1) x))
          
          ;; Establish starting conditions for derivatives
          ;; (NOTE: default case already has u and y2 zeroed.
          (case kind
            (:natural)
            (:clamped
             (setf (aref y2 0) -0.5)
             (if yder1
                 (let* ((dx (- (aref x 1) (aref x 0)))
                        (dy (- (aref y 1) (aref y 0))))
                   (declare (type single-float dx dy))
                   (setf (aref u 0) (* (/ 3.0 dx)
                                       (- (/ dy dx) yder1))
                         )))))
          
          ;; Invert the tridiagonal system of equations to solve
          ;; for node second derivatives
          (loop for i from 1 below nlast do
                (let* ((dx10 (- (aref x i) (aref x (1- i))))
                       (dy10 (- (aref y i) (aref y (1- i))))
                       (dx21 (- (aref x (1+ i)) (aref x i)))
                       (dy21 (- (aref y (1+ i)) (aref y i)))
                       (dx20 (+ dx21 dx10))
                       (sig  (/ dx10 dx20))
                       (p    (+ 2.0 (* sig (aref y2 (1- i)))))
                       (ui   (- (/ dy21 dx21)
                                (/ dy10 dx10))))
                  
                  (declare (type fixnum i)
                           (type single-float dx10 dy10 dx21 dy21 dx20 sig p ui))
                  
                  (setf (aref y2 i) (/ (- sig 1.0) p)
                        (aref u i)  (/ (- (/ (* 6.0 ui) dx20)
                                          (* sig (aref u (1- i))))
                                       p)
                        )))
        
          ;; establish the final boundary conditions on the equations
          (case kind
            (:natural)
            (:clamped
             (let* ((qn 0.5)
                    (dx (- (aref x nlast) (aref x (1- nlast))))
                    (un (if yder2
                            (let ((dy (- (aref y nlast) (aref y (1- nlast)))))
                              (* (/ 3.0 dx)
                                 (- yder2 (/ dy dx))))
                          0.0)))
               
               (declare (type single-float qn dx un))
               
               (setf (aref y2 nlast) (/ (- un (* qn (aref u nlast)))
                                        (+ 1.0 (* qn (aref y2 nlast))))
                     ))))
          
          ;; backsubstitution for solution of node second derivatives
          (loop for i from (1- nlast) downto 0 do
                (locally
                  (declare (type fixnum i))
                  (setf (aref y2 i) (+ (* (aref y2 i) (aref y2 (1+ i)))
                                       (aref u i)))))
          
          ;; return an interpolator suitable for use in SPLINT
          (make-instance '<spline-interpolator>
                         :xtable  x
                         :ytable  y
                         :y2table y2
                         :nlast   nlast)
          )))))
  
(defmethod splint ((x integer) (interp-record <spline-interpolator>))
  ;; ensure lookup value is a float
  (splint (sfloat x) interp-record))

(defmethod splint ((x double-float) (interp-record <spline-interpolator>))
  (splint (sfloat x) interp-record))

(defmethod splint ((x single-float) (interp-record <spline-interpolator>))
  ;; Cubic-spline interpolation using a precomputed <SPLINE-INTERPOLATOR>
  ;; as obtained from an initial call to routine SPLINE.
  ;;
  ;; If the lookup value X lies beyond the domain of the table,
  ;; then simple linear extrapolation is provided.
  ;; [NOTE: This means that if you expect to have X values outside of the table
  ;;  domain, AND if you expect continuous first dirivatives in Y, then you should
  ;;  construct the <SPLINE-INTERPOLATOR> using keyword :KIND = :EXTENSIONAL].
  ;;
  ;; Adapted from "Numerical Recipes in C - Second Edition",
  ;; by Press, Teukolsky, Vetterling, and Flannery,
  ;; Cambridge University Press, 1996, Section 3.3 and 3.4, pp. 113-117.
  ;;
  ;; DM/MCFA  04/02
  ;; -------------------------------------------------
  ;; cubic spline interpolation using a previously computed
  ;; <spline-interpolator> object as obtained from an initial
  ;; call to SPLINE.

  (with-accessors ((spl-x0-last  spl-x0-last)
                   (spl-x1-last  spl-x1-last)
                   (spl-y0-last  spl-y0-last)
                   (spl-y1-last  spl-y1-last)
                   (spl-y20-last spl-y20-last)
                   (spl-y21-last spl-y21-last)
                   (spl-h-last   spl-h-last)
                   (spl-xtable   spl-xtable)
                   (spl-ytable   spl-ytable)
                   (spl-y2table  spl-y2table)
                   (spl-nlast    spl-nlast)) interp-record
    
    ;; if x between previous lookup nodes then quick interpolation
    ;; with already located node data.
    (if (and spl-x0-last
             spl-x1-last
             (<= spl-x0-last x spl-x1-last))
        (let* ((h spl-h-last)
               (a (/ (- spl-x1-last x) h))
               (b (/ (- x spl-x0-last) h)))
          
          (declare (type single-float h a b))
          
          (+ (* a spl-y0-last)
             (* b spl-y1-last)
             (/ (* (+ (* (- (* a a a) a) spl-y20-last)
                      (* (- (* b b b) b) spl-y21-last))
                   h h)
                6.0))
          ))
    
    ;; else -- binary search for bracketing nodes
    (let* ((xtable  spl-xtable)  ;; dereference from object for speedup
           (ytable  spl-ytable)
           (y2table spl-y2table)
           (nlast   spl-nlast))
      
      (declare (type (array single-float (*)) xtable ytable y2table)
               (type fixnum nlast))
      
      (cond
       ((< x (aref xtable 0))
        ;; off the fore end of the table -- use linear extrapolation
        (let* ((h (- (aref xtable 1) (aref xtable 0)))
               (a (/ (- (aref xtable 1) x) h))
               (b (/ (- x (aref xtable 0)) h)))
          (declare (type single-float h a b))
          (+ (* a (aref ytable 0))
             (* b (aref ytable 1)))
          ))
       
       ((> x (aref xtable nlast))
        ;; off the aft end of the table -- use linear extrapolation
        (let* ((h (- (aref xtable nlast) (aref xtable (1- nlast))))
               (a (/ (- (aref xtable nlast) x) h))
               (b (/ (- x (aref xtable (1- nlast))) h)))
          (declare (type single-float h a b))
          (+ (* a (aref ytable (1- nlast)))
             (* b (aref ytable nlast)))
          ))
       
       (t  ;; inside the table limits -- use cubic spline interpolation
           ;; locate surrounding nodes with a binary search
           (let ((n0  0)
                 (n1  nlast))
             (declare (type fixnum n0 n1))
             (loop while (> (- n1 n0) 1) do
                   (let* ((nmid (truncate (+ n0 n1) 2)))
                     (declare (type fixnum nmid))
                     (if (< x (aref xtable nmid))
                         (setf n1 nmid)
                       (setf n0 nmid))
                     ))
             
             ;; form the cubic spline interpolate and also cache node info
             ;; in case next lookup is incremental and inside the same pair
             ;; of surrounding nodes.
             (let* ((x0  (aref xtable n0))
                    (x1  (aref xtable n1))
                    (y0  (aref ytable n0))
                    (y1  (aref ytable n1))
                    (y20 (aref y2table n0))
                    (y21 (aref y2table n1))
                    (h  (- x1 x0))
                    (a  (/ (- x1 x) h))
                    (b  (/ (- x x0) h)))
               
               (declare (type single-float x0 x1 y0 y1 y20 y21 h a b))
               
               ;; cache for possible speedup on next lookup
               (setf spl-x0-last  x0
                     spl-x1-last  x1
                     spl-y0-last  y0
                     spl-y1-last  y1
                     spl-y20-last y20
                     spl-y21-last y21
                     spl-h-last   h)
               
               ;; compute return interoplated value
               (+ (* a y0)
                  (* b y1)
                  (/ (* (+ (* (- (* a a a) a) y20)
                           (* (- (* b b b) b) y21))
                        h h)
                     6.0))
               )))
       ))))

                          
#|
(defparameter xtable
  #(1 2.4 3.3 5.5 7.9 10))
(defparameter ytable
  #(2 4.8 6.6 11 15.8 20))

(let* ((ytable (reverse ytable))
       (x  (vmath:vectorwise ((x (vm:framp 150)))
                             (/ x 10)))
       (spl (spline xtable ytable
                    :kind :extensional))
       (y3 (vmath:vectorwise ((x x))
                             (splint x spl)))
       (y1 (vmath:vectorwise ((x x))
                            (interp x xtable ytable
                                    :order 1)))
       (y2 (vmath:vectorwise ((x x))
                            (interp x xtable ytable
                                    :order 2))))
  (sg:plot x y1 :color sg:$red
           :xrange '(0 2)
           :yrange '(15 25)
           )
  (sg:oplot x y2 :color sg:$darkgreen)
  (sg:oplot x y3 :color sg:$blue)
  (sg:oplot xtable ytable :symbol sg:$sym-circle)
  )

;; This version gives completely different quadratic interpolation results
;; since the intervals are reversed. But linear and cubic-spline interpolation
;; produces the same results as before...
(let* ((xtable (reverse xtable))
       (x  (vmath:vectorwise ((x (vm:framp 150)))
                             (/ x 10)))
       (y1 (vmath:vectorwise ((x x))
                            (interp x xtable ytable
                                    :order 1)))
       (y2 (vmath:vectorwise ((x x))
                            (interp x xtable ytable
                                    :order 2)))
       (spl (spline xtable ytable
                    :kind :extensional))
       (y3 (vmath:vectorwise ((x x))
                             (splint x spl))))
  (sg:plot x y1 :color sg:$red
           :xrange '(0 2)
           :yrange '(15 25)
           )
  (sg:oplot x y2 :color sg:$darkgreen)
  (sg:oplot x y3 :color sg:$blue)
  (sg:oplot xtable ytable :symbol sg:$sym-circle)
  )

;; This produces a zero result (as it should) for the difference between
;; cubic spline interpolation between the two cases of reversed X and reversed Y.
(let* ((x  (vmath:vectorwise ((x (vm:framp 150)))
                             (/ x 10)))
       (spl1 (spline (reverse xtable) ytable
                    :kind :extensional))
       (spl2 (spline xtable (reverse ytable)
                    :kind :extensional))
       (y1 (vmath:vectorwise ((x x))
                             (splint x spl1)))
       (y2 (vmath:vectorwise ((x x))
                             (splint x spl2)))
       (dy (vmath:vectorwise ((y1 y1)
                              (y2 y2))
                             (- y1 y2))))

  (sg:wset 2)
  (sg:plot x dy)
  (sg:wset 1)
  (sg:plot x y1 :color sg:$red
            :xrange '(0 5)
            :yrange '(10 25)
            )
  (sg:oplot x y2 :color sg:$darkgreen)
  (sg:oplot (reverse xtable) ytable :symbol sg:$sym-circle)
  (sg:oplot xtable (reverse ytable) :symbol sg:$sym-box)
  )
|#

(defun interp (x xtable ytable &key (order 1))
  ;; Linear and quadratic interpolation of tables.
  ;; The XTABLE can be in either ascending or descending order.
  ;; It makes no difference for linear interpolation,
  ;; but it does make a slight difference for quadratic interpolation.
  ;;
  ;; Unless you just want a simple linear interpolation, it is recommended
  ;; to use SPLINT for cubic-spline interpolation. That routine presents
  ;; results with continuous first derivatives. This routine cannot.
  ;;
  ;; DM/MCFA 04/02
  ;; ---------------------------------------------------------------
  (let* ((nel   (length xtable))
         (nlast (1- nel))
         (n0    0)
         (n1    nlast)
         (x0    (aref xtable n0))
         (x1    (aref xtable n1))
         (le    (if (> x1 x0)
                    #'<=
                  #'>=)))
    (assert (> nel order))
    (assert (= nel (length ytable)))
    (cond
     ((= x x0) (aref ytable n0))
     ((= x x1) (aref ytable n1))
     ((funcall le x0 x x1)
      (let* ((nmid (truncate (+ n0 n1) 2))
             (xmid (aref xtable nmid)))
        (do ()
            ((or (= x xmid)
                 (= n0 nmid)))
          (if (funcall le x xmid)
              (setf n1 nmid
                    x1 xmid)
            (setf n0 nmid
                  x0 xmid))
          (setf nmid (truncate (+ n0 n1) 2)
                xmid (aref xtable nmid)))
        (if (= x xmid)
            (aref ytable nmid)
          (progn
            (assert (= (1+ n0) n1))
            (ecase order
              (1
               (let ((frac (/ (- x x0) (- x1 x0))))
                 (+ (* frac (aref ytable n1))
                    (* (- 1 frac) (aref ytable n0)))
                 ))
              (2
               (let* ((dx0 (- x x0))
                      (dx1 (- x x1))
                      (y0  (aref ytable n0))
                      (y1  (aref ytable n1))
                      (n2  (cond
                            #|
                            ;; Choose 3rd node closest to datum.
                            ;; NOT RECOMMENDED -- it produces abrupt
                            ;; discontinuities at the points halfway between
                            ;; fitting nodes...
                            ((zerop n0) 2)
                            ((= n1 nlast) (- nlast 2))
                            (t (if (<= (abs dx0) (abs dx1))
                                   (1- n0)
                                 (1+ n1)))
                            |#
                            ;; Always choose 3rd node on the high side
                            ;; except for the very first interval.
                            ((= n1 nlast) (- nlast 2))
                            (t (+ n0 2))
                            ))
                      (x2  (aref xtable n2))
                      (dx2 (- x x2))
                      (y2  (aref ytable n2))
                      (dx10 (- x1 x0))
                      (dx20 (- x2 x0))
                      (dx21 (- x2 x1))
                      (dx12 (- dx21)))
                 (+ (/ (* y0 dx1 dx2) (* dx10 dx20))
                    (/ (* y1 dx0 dx2) (* dx10 dx12))
                    (/ (* y2 dx0 dx1) (* dx20 dx21)))
                 ))
              )))
        ))
     ((funcall le x1 x)
      ;; extrapolation from one end -- always linear
      (let* ((x0   (aref xtable (1- n1)))
             (frac (/ (- x x0)
                      (- x1 x0))))
        (+ (* frac (aref ytable n1))
           (* (- 1 frac) (aref ytable (1- n1))))
        ))
     (t
      ;; extrapolation from the other end -- always linear
      (let* ((frac (/ (- x x0)
                      (- (aref xtable (1+ n0)) x0))))
        (+ (* frac (aref ytable (1+ n0)))
           (* (- 1 frac) (aref ytable n0)))
        ))
     )))

#|
(defparameter xtable
  #(1 2.4 3.3 5.5 7.9 10))
(defparameter ytable
  #(2 4.8 6.6 11 15.8 20))


(interp 0.5 (reverse xtable) ytable)
(interp 1.1 #(1 2) #(3 4)) ;; test smallest valid tables
(interp 1.1 #(1) #(3))  ;; should generate assertion violation

(let* (;(xtable #(1 2 3 4 5 6))
       ;(ytable #(2 4 6 8 10 12))
       (x  (vmath:vectorwise ((x (vm:framp 150)))
                             (/ x 10)))
       (y1 (vmath:vectorwise ((x x))
                            (interp x (reverse xtable) ytable
                                    :order 1)))
       (y2 (vmath:vectorwise ((x x))
                            (interp x (reverse xtable) ytable
                                    :order 2))))
  #|
  (sg:plot x y1
           :symbol sg:$sym-circle)
  (sg:oplot x y2
            :symbol sg:$sym-box)
  |#
  (sg:plot x y1 :color sg:$red)
  (sg:oplot x y2 :color sg:$darkgreen)
  (sg:oplot (reverse xtable) ytable :symbol sg:$sym-circle)
  )

|#
