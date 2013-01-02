;; Tchebyshev.lisp -- Various useful Tchebyshev Interpolators/Approximators
;; DM/RAL  12/06

(defclass <approximation> ()
  ((xmin                     :accessor approximation-xmin
                             :initarg  :xmin)
   (xmax                     :accessor approximation-xmax
                             :initarg  :xmax)
   (range-reduction-function :accessor approximation-range-reduction-function
                             :initarg :range-reduction-function)
   (coefficients             :accessor approximation-coefficients
                             :initarg :coefficients)))

(defclass <power-series-approximation>    (<approximation>) ())
(defclass <tchebyshev-approximation>      (<approximation>) ())
(defclass <tchebyshev-even-approximation> (<tchebyshev-approximation>) ())

;; ---------------------------------------------------------------------------
;; functional folding operators
;; Thesse fold over the integers from 0 to N-1, in forward (foldil)
;; and reverse (foldir) order.
;;
;; (FoldIR N Fn Init) folds function Fn from index N-1 downto 0,
;; with Init as the initial state. Fn takes args Index and State and returns a new state.
(defun foldir (n fn init)
  (um:with-tail-pure-code
    (labels ((iter (ix state)
               (if (minusp ix)
                   state
                 (iter (1- ix) (funcall fn ix state)))))
      (iter (1- n) init))))

;; (FoldIL N Fn Init) folds function Fn from index 0 to N-1,
;; with Init as the initial state. Fn takes args State and Index and returns a new state.
(defun foldil (n fn init)
  (um:with-tail-pure-code
    (labels ((iter (state ix)
               (if (>= ix n)
                   state
                 (iter (funcall fn state ix) (1+ ix)))))
      (iter init 0))))

;; E.g., (foldil (lambda (ix sum) (+ ix sum)) 10) ==> 45

;; ---------------------------------------------------------------------------
(defmethod make-approximator ((approx <power-series-approximation>))
  (with-accessors ((rrfn  approximation-range-reduction-function)
                   (coffs approximation-coefficients)) approx
    (let* ((order (1- (length coffs))))
      #'(lambda (x)
          (let ((u (funcall rrfn x)))
            (foldir order
                    (lambda (ix y)
                      (+ (* u y) (aref coffs ix)))
                    (aref coffs order))
            ))
      )))

(defmethod make-approximator ((approx <tchebyshev-approximation>))
  (with-accessors ((rrfn  approximation-range-reduction-function)
                   (coffs approximation-coefficients)) approx
    (let* ((order (1- (length coffs))))
      #'(lambda (x)
          (let ((u (funcall rrfn x)))
            (destructuring-bind (p1 p2)
                (foldir order
                        (lambda (ix state)
                          (destructuring-bind (p1 p2) state
                            (list (+ (- (* 2.0d0 u p1) p2) (aref coffs ix)) p1)))
                        (list (aref coffs order) 0.0d0))
              (- p1 (* u p2) (* 0.5d0 (aref coffs 0)))
              )))
      )))

;; ---------------------------------------------------------------------------
(defun inverse-range-reduction-function (xmin xmax)
  (lambda (x)
    (* 0.5d0 (+ (* x (- xmax xmin)) xmax xmin))))

(defun range-reduction-function (xmin xmax)
  (lambda (x)
    (/ (- (* 2.0d0 x) xmin xmax) (- xmax xmin))))


;; produce a Tchebyshev approximation in T(n,x) for x in [-1,1]
(defun make-tchebyshev-interpolation (&key order xmin xmax function)
  (let* ((zs    (loop for ix from 0 to order collect
                      (/ (* pi (1+ (* 2 ix))) (* 2 (1+ order)))))
         (xs    (mapcar (inverse-range-reduction-function xmin xmax)
                        (mapcar #'cos zs)))
         (ys    (mapcar function xs))
         (coffs (make-array (1+ order) :initial-element 0.0d0))
         (sf    (/ 2.0d0 (1+ order))))
    (loop for ix from 0 to order do
          (setf (aref coffs ix)
                (* sf (loop for y in ys
                            and z in zs
                            sum
                            (* y (cos (* ix z)))))
                ))
        (make-instance '<tchebyshev-approximation>
                       :xmin                     xmin
                       :xmax                     xmax
                       :range-reduction-function (range-reduction-function xmin xmax)
                       :coefficients             coffs)
        ))

;; produce a power-series approximation in x for x in [-1,1]
(defmethod convert-tchebyshev-to-power-series-approximation
           ((approx <tchebyshev-approximation>))
  (let* ((coffs (approximation-coefficients approx))
         (ord   (1- (length coffs))))
    (labels ((make-zeros ()
               (make-array (1+ ord) :initial-element 0.0d0)))
      (destructuring-bind (p1 p2)
          (foldir (1+ ord)
                  (lambda (ix state)
                    (destructuring-bind (p1 p2) state
                      (let* ((p0 (make-zeros)))
                        ;; p0 = 2*x*p1 - p2 + c[n]
                        ;; the multiply by x is equivalent to a right shift of the p1 array
                        ;;
                        ;; So for index 0: c[n] - p2[0]
                        ;; and for others: 2 p1[ix-1] - p2[ix]
                        (setf (aref p0 0) (- (aref coffs ix)
                                             (aref p2 0)))
                        (loop for jx from 1 to ord do
                              (setf (aref p0 jx)
                                    (- (* 2.0d0 (aref p1 (1- jx)))
                                       (aref p2 jx))))
                        (list p0 p1))))
                  (list (make-zeros) (make-zeros)))
        (let ((p0 (make-zeros)))
          ;; final step: ans = p1 - x*p2 - c[0]/2
          ;; the multiply by x is equivalent to a right shift of the p1 array
          ;;
          ;; So for index 0: p1[0] - c[0]/2
          ;; and for others: p1[ix] - p2[ix-1]
          (setf (aref p0 0)
                (- (aref p1 0)
                   (* 0.5d0 (aref coffs 0))))
          (loop for ix from 1 to ord do
                (setf (aref p0 ix)
                      (- (aref p1 ix)
                         (aref p2 (1- ix)))))
          (make-instance '<power-series-approximation>
                         :xmin
                         (approximation-xmin approx)
                         :xmax
                         (approximation-xmax approx)
                         :range-reduction-function
                         (approximation-range-reduction-function approx)
                         :coefficients p0)
          )))))

;; ---------------------------------------------------------------------------
#|              
;; Example: approximation of Exp(x) for x in [-1,1]
;; Actually, we approximate the function P(x) = (Exp(x) - 1)/x so that our
;; approximation will have zero error at x = 0, exp(x) = 1 + x P(x)
(let* ((order  9)
       (xmin -1.0d0)
       (xmax  1.0d0)
       (apprx  (make-tchebyshev-interpolation
                :order order
                :xmin xmin
                :xmax xmax
                :function (lambda (x)
                            (if (zerop x)
                                1.0d0
                              (/ (- (exp x) 1.0d0) x)))))
       (papprx (convert-tchebyshev-to-power-series-approximation apprx))
       (dom (list xmin xmax))
       (approx (make-approximator  apprx))
       (papprox (make-approximator papprx))
       (err (lambda (x)
              (- (+ 1.0d0 (* x (funcall approx x))) (exp x))))
       (errp (lambda (x)
               (- (+ 1.0d0 (* x (funcall papprox x))) (exp x)))))
  ;; (inspect apprx)
  ;; (inspect papprx)
  (let ((win (plt:wset 0 :clear t)))
    (plt:fplot win dom err)
    (plt:fplot win dom errp :color :red)
    
    ;; show that the unconstrained approximation has lower peak errors
    ;; BTW, the unconstrained approximation is already producing zero error at x = 0,
    ;; as long as its order is even, so constraining in this case is pointless.
    (let* ((order  (+ order 1))
           (apprx (make-tchebyshev-interpolation
                   :order order
                   :xmin  xmin
                   :xmax  xmax
                   :function #'exp))
           (approx (make-approximator apprx))
           (papprx (convert-tchebyshev-to-power-series-approximation apprx))
           (err (lambda (x)
                  (- (funcall approx x) (exp x)))))
      ;; (inspect papprx)
      (plt:fplot win dom err :color :blue))))

(let* ((order 2)
       (xmin -1.0d0)
       (xmax  1.0d0)
       (apprx (make-tchebyshev-interpolation
               :order order
               :xmin  xmin
               :xmax  xmax
               :function #'exp))
       (approx (make-approximator apprx))
       (dom  (list xmin xmax))
       (err (lambda (x)
              (- (funcall approx x) (exp x)))))
  (let ((win (plt:wset 0 :clear t)))
    (plt:fplot win dom err :color :blue)))


;; Example: approximation of Log(2,x) for x in [0.5,1]
;; Actually approximate log(2,x) = (1-x)*(-2 + (x-0.5)*P(x))
;; where P(x) = (Log(2,x)+2(1-x))/((1-x)(x-0.5))
;; so that we have zero error at x = 1 and at x = 0.5.
(let* ((order   6)
       (approx (make-tchebyshev-interpolation
                :order order
                :xmin  0.5
                :xmax  1
                :function (lambda (x)
                            (cond
                             ((= x 1.0d0) (- 4.0d0 (/ 2.0d0 (log 2.0d0))))
                             ((= x 0.5d0) (- (/ 4.0d0 (log 2.0d0)) 4.0d0))
                             (t (/ (+ (log x) (* 2.0d0 (- 1.0d0 x) (log 2.0d0)))
                                   (* (log 2.0d0) (- 1.0d0 x) (- x 0.5d0))))
                             ))
                ))
       (papprox (convert-tchebyshev-to-power-series-approximation approx))
       (x '(0.5d0 1.0d0))
       (approx (make-approximator  approx))
       (papprox (make-approximator papprox))
       (err (lambda (x)
              (- (* (- 1.0d0 x)
                    (- (* (- x 0.5d0) (funcall approx x))
                       2.0d0))
                 (log x 2.0d0))))
       (errp (lambda (x)
               (- (* (- 1.0d0 x)
                     (- (* (- x 0.5d0) (funcall papprox x))
                        2.0d0))
                  (log x 2.0d0)))))
  (let ((win (plt:wset 1 :clear t)))
    (plt:fplot win x err)
    (plt:fplot win x errp :color :red))
  
  (let ((win (plt:wset 2 :clear t)))
    (plt:fplot win x approx))
  
  (let* ((order  (+ order 2))
         (apprx (make-tchebyshev-interpolation
                 :order order
                 :xmin  0.5
                 :xmax  1
                 :function (lambda (x) (log x 2.0d0))))
         (approx (make-approximator apprx))
         (err (lambda (x)
                (- (funcall approx x) (log x 2.0d0)))))
    (let ((win (plt:wset 1 :clear t)))
      (plt:fplot win x err :color :blue))))

|#

;; ---------------------------------------------------------------------------
;; Special versions for Even Functions where Fn(-x) = Fn(x)
;;
(defun even-range-reduction-function (xmin xmax)
  (lambda (x)
    (/ (- x xmin) (- xmax xmin))))

(defun inverse-even-range-reduction-function (xmin xmax)
  (lambda (x)
    (+ (* x (- xmax xmin)) xmin)))

(defun sqr (x)
  (* x x))

(defun twoxm1 (x)
  (- (* 2.0d0 x) 1.0d0))

;; For even functions fn(-x) = fn(x) and all odd T(2n+1,x) = 0.
;; Produce a Tchebyshev approximation in T(2n,x) = T*(n,x^2) for x in [0,1]
;; based on definition of shifted Tchebyshev Polynomials T*(n,x) = T(n,2x-1)
;; Tchebyshev polynomials T(n,x) are equal-ripple over domain x in [-1,1]
;; Shifted Tchebyshev polynomials T*(n,x) are equal-ripple over domain x in [0,1]
;; where:  T(n,T(m,x)) = T(n*m,x)
;; and:    T*(x^2) = T(2n,x)
;;
(defun make-even-tchebyshev-interpolation (&key order xmin xmax function)
  (let* ((zs    (loop for ix from 0 to (* 2 order) collect
                      (/ (* pi (1+ (* 2 ix))) (* 2 (1+ (* 2 order))))))
         (xs    (mapcar (inverse-even-range-reduction-function xmin xmax)
                        (mapcar #'abs (mapcar #'cos zs))))
         (ys    (mapcar function xs))
         (coffs (make-array (1+ order) :initial-element 0.0d0))
         (sf    (/ 2.0d0 (1+ (* 2 order)))))
    (loop for ix from 0 to order do
          (setf (aref coffs ix)
                (* sf (loop for y in ys
                            and z in zs
                            sum
                            (* y (cos (* 2 ix z)))))
                ))
        (make-instance '<tchebyshev-even-approximation>
                       :xmin                     xmin
                       :xmax                     xmax
                       :range-reduction-function (um:compose
                                                  #'twoxm1 #'sqr
                                                  (even-range-reduction-function xmin xmax))
                       :coefficients             coffs)
        ))

;; produce a power series approximation in x^2 for x in [0,1]
(defmethod convert-tchebyshev-to-power-series-approximation
           ((approx <tchebyshev-even-approximation>))
  (let* ((coffs (approximation-coefficients approx))
         (ord   (1- (length coffs))))
    (labels ((make-zeros ()
               (make-array (1+ ord) :initial-element 0.0d0)))
      (destructuring-bind (p1 p2)
          (foldir (1+ ord)
                  (lambda (ix state)
                    (destructuring-bind (p1 p2) state
                      (let* ((p0 (make-zeros)))
                        ;; p0 = 2*(2*x-1)*p1 - p2 + c[n]
                        ;; in other words... same as for regular Tchebyshev approx
                        ;; but where x is replaced by (2x-1)
                        ;;
                        ;; So, for index 0: c[n] - 2 p1[0] - p2[0]
                        ;; and for others:  4 p1[ix-1] - 2 p1[ix] - p2[ix]
                        (setf (aref p0 0) (- (aref coffs ix)
                                             (* 2.0d0 (aref p1 0))
                                             (aref p2 0)))
                        (loop for jx from 1 to ord do
                              (setf (aref p0 jx)
                                    (- (* 4.0d0 (aref p1 (1- jx)))
                                       (* 2.0d0 (aref p1 jx))
                                       (aref p2 jx))))
                        (list p0 p1))))
                  (list (make-zeros) (make-zeros)))
        (let ((p0 (make-zeros)))
          ;; final step: ans = p1 - (2*x-1)*p2 - c[0]/2
          ;; or, for index 0:  p1[0] + p2[0] - c[0]/2
          ;; and for others:   p1[ix] + p2[ix] - 2 p2[ix-1]
          (setf (aref p0 0)
                (- (+ (aref p1 0)
                      (aref p2 0))
                   (* 0.5d0 (aref coffs 0))))
          (loop for ix from 1 to ord do
                (setf (aref p0 ix)
                      (- (+ (aref p1 ix)
                            (aref p2 ix))
                         (* 2.0d0 (aref p2 (1- ix))))
                      ))
          (make-instance '<power-series-approximation>
                         :xmin
                         (approximation-xmin approx)
                         :xmax
                         (approximation-xmax approx)
                         :range-reduction-function
                         (um:compose #'sqr
                                  (even-range-reduction-function
                                   (approximation-xmin approx)
                                   (approximation-xmax approx)))
                         :coefficients p0)
          )))))


;; ---------------------------------------------------------------------------
#|
(let* ((apprx (make-even-tchebyshev-interpolation
                :order 4
                :xmin  0
                :xmax  1
                :function (lambda (x)
                            (if (zerop x) (* 0.5d0 pi)
                              (/ (sin (* 0.5d0 pi x)) x)))))
       (papprx (convert-tchebyshev-to-power-series-approximation apprx))
       (x '(0.0d0 1.0d0))
       (approx (make-approximator  apprx))
       (papprox (make-approximator papprx))
       (err (lambda (x)
              (- (* x (funcall approx x)) (sin (* 0.5d0 pi x)))))
       (errp (lambda (x)
               (- (* x (funcall papprox x)) (sin (* 0.5d0 pi x))))))
  ;; (inspect apprx)
  ;; (inspect papprx)
  (let ((win (plt:wset 0 :clear t)))
    (plt:fplot win x err)
    (plt:fplot win x errp :color :red))
  (inspect (approximation-coefficients papprx)))
|#
;; ---------------------------------------------------------------------------

#|
(let* ((delta 0.05d0)
       (fn    (lambda (x)
                (let ((x (tan (* pi x 0.5d0))))
                  (+ (* x x) (log (erfc x))))))
       (apprx (make-tchebyshev-interpolation
               :order 11
               :xmin  0.0d0
               :xmax  (- 1.0d0 delta)
               :function fn))
       (papprx (convert-tchebyshev-to-power-series-approximation apprx))
       (approx (make-approximator apprx))
       (papprox (make-approximator papprx))
       (dom   (list 0.0d0 (- 1.0d0 delta)))
       (err   (lambda (x)
                (- (funcall approx x) (funcall fn x))))
       (perr  (lambda (x)
                (- (funcall papprox x) (funcall fn x)))))
  ;(inspect apprx)
  (let ((win (plt:wset 1 :clear t)))
    (plt:fplot win dom err)
    ;;(plt:fplot win dom perr :color :red)
    )
  )




;; -----------------------------------------
             
|#