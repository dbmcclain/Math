;; divulge.lisp -- William Kahan's checkout routines for Floating Point
;; performance of machines and languages
;;
;; DM/MCFA  06/04

(proclaim '(debug 2)) ;; eliminate tail recursion

#-:WIN32
(fli:register-module "/usr/lib/libm.dylib")

#-:WIN32
(fli:define-foreign-function (nextafter "nextafter" :source)
    ((x :double)
     (y :double))
  :result-type :double
  :language :ansi-c)

#+:WIN32
(fli:define-foreign-function (nextafter "_nextafter" :source)
    ((x :double)
     (y :double))
  :result-type :double
  :language :ansi-c)

#-:WIN32
(fli:define-foreign-function (nextafterf "nextafterf" :source)
    ((x :float)
     (y :float))
  :result-type :float
  :language :ansi-c)

(defparameter eps (- (nextafter 1.0D0 2.0D0) 1.0D0))

(defun pnl (item)
  (princ item)
  (terpri))

(defun v<*> (v1 v2)
  (cond
   ((consp v1)
    (cond
     ((consp v2)
      (labels ((vmul (v1 v2 ans)
                 (cond
                  ((and v1 v2)
                   (vmul (rest v1)
                         (rest v2)
                         (+ ans (* (first v1) (first v2)))
                        ))
                  ((or v1 v2)
                   (error "v<*>: vectors must have same length"))
                  (t ans))
                 ))
        (vmul v1 v2 0)))
     (t
      (mapcar #'(lambda (v)
                  (* v v2))
              v1))
     ))
   ((consp v2)
    (mapcar #'(lambda (v)
                (* v v1))
            v2))
   (t (* v1 v2))
   ))

(defun reduce-rank (mat)
  (cond
   ((consp mat)
    (cond
     ((= 1 (length mat))
      (reduce-rank (first mat)))
     (t mat)))
   (t mat)))

(defun transpose (mat)
  (cond
   ((consp mat)
    (cond
     ((consp (first mat))
      (labels ((trn (mat ans)
                 (cond
                  ((consp (first mat))
                   (trn (mapcar #'rest mat)
                        (cons (mapcar #'first mat) ans)))
                  (t (nreverse ans))
                  )))
        (reduce-rank (trn mat nil))
        ))
     (t (mapcar #'list mat))
     ))
   (t mat)
   ))
    
(defun m* (m1 m2)
  (labels ((mrcs (r cs ans)
             (cond
              ((consp (first cs))
               (mrcs r (rest cs) (cons (v<*> r (first cs)) ans)))
              (cs  (v<*> r cs))
              (t  (reduce-rank (nreverse ans)))))
           (mrscs (rs cs ans)
             (cond
              ((consp (first rs))
               (mrscs (rest rs) cs (cons (mrcs (first rs) cs nil) ans)))
              (rs (mrcs rs cs nil))
              (t  (reduce-rank (nreverse ans))))))
    (mrscs m1 (transpose m2) nil)
    ))

(defun vmax (v)
  (labels ((vm (v ans)
             (if v
                 (vm (rest v) (max (first v) ans))
               ans)))
    (vm v 0)))

(defun t1 ()
  (pnl (format nil "~A ~A version ~A"
               (machine-type)
               (lisp-implementation-type)
               (lisp-implementation-version)))
  (pnl "Machine: eMac 1.25 GHz G4"))

(defun t2 ()
  (pnl (if (zerop (v<*> `(,eps 9.0d0 9.0d0)
                        `(,eps 9.0d0 -9.0d0)))
           "Scalar product is summed from left to right"
         "Scalar product is summed from right to left")))

(defun t3 ()
  (let* ((mxmuleps (let* ((z (/ 4.0d0 3.0d0))
                          (u (+ (* 3.0d0 (- 1.0d0 z)) 1.0d0)))
                     (unless (equalp (abs u) eps)
                       (pnl "Has precision been altered?")
                       (pnl "Why do the following differ...")
                       (pnl (format nil "anULPofOne = ~F, Eps = ~F"
                                    (abs u) eps))
                       (pnl "Now mxmuleps cannot be trusted"))
                     (let* ((z    (* 4.0d0 (- z 1.0d0)))
                            (zzz `((,z) (,z) (,z)))
                            (mat `((,(- u) ,(+ 1.0d0 u) -1.0d0)
                                   (-1.0d0 ,(+ 1.0d0 u) ,(- u))))
                            (y (* 3.0d0 (vmax
                                     (mapcar #'abs 
                                             (m* mat zzz))
                                     ))))
                       (if (zerop y)
                           'Nan
                         (/ u (round (/ u y))))
                       )))
         (y (v<*> `(,(- 1.0d0 eps) ,(- eps 1.0d0)) 
                  `(,(+ 1.0d0 eps) ,(+ eps 1.0d0))))
         (x (v<*> `(,eps 9.0d0 9.0d0) 
                  `(,eps 9.0d0 -9.0d0))))
    (multiple-value-bind (me yprime)
        (case mxmuleps
          (NaN (values " is NaN "
                       (or (zerop y)
                           (not (equalp x
                                        (if (minusp y) 1 0))))))
          (otherwise (if (equalp mxmuleps eps)
                         (values " = eps " y)
                       (values (format nil " - eps/~D"
                                       (/ eps mxmuleps))
                               y)))
          )
      (declare (ignore yprime))
      (pnl (format nil " mxmuleps~A" me))
      (unless (zerop y)
        (pnl "Something about the previous two statements is wrong."))
      )))
               
         
(defun t4 ()
  (let* ((x (loop for ix from 1 to 128 collect
                  (+ 1.0d0 (/ (float ix 1.0d0) 65536.0d0))))
         (n (loop for ix from 1 to 128 collect
                  (+ 65500.0d0 (* 2.0d0 (float ix 1.0d0)))))
         (x3n (loop for x in x and n in n collect
                    (expt x (* 3.0d0 n))))
         (xxxn (loop for x in x and n in n collect
                     (expt (* x x x) n)))
         (diff (/ (vmax (mapcar #'(lambda (x3n xxxn)
                                    (abs (/ (- x3n xxxn) x3n)))
                                x3n xxxn)) 
                  eps)))
    (pnl (format nil "y^N errors as large as ~F ULPs have been seen."
                 diff))
    ))

(defun unoise (&optional nel)
  (if nel
      (loop for ix from 0 below nel collect
            (random 1.0d0))
    (random 1.0d0)))

(defun t5 ()
  (let* ((x (loop for ix from 0 below (* 32 32) collect
                  (+ 0.5d0 (unoise))))
         (n (loop for ix from 0 below (* 32 32) 
                  and x in x collect
                  (+ 7.0d0 (round (* 256.0d0 x)))))
         (y2xx (vmax (mapcar #'(lambda (x n)
                                 (- (abs (expt x (* 2.0d0 n)))
                                    (expt (* x x) n)))
                             x n)))
         (notso (if (zerop y2xx) "not " "")))
    (princ "y^N is ")
    (princ notso)
    (princ "accumulated extra-precisely.")
    (terpri)))

(defun t6 ()
  (let* ((y (+ 1.0d0 eps))
         (d (- 1.0d0 (/ eps 2.0d0)))
         (x (+ 1.0d0 (* y (/ eps 2.0d0))))
         (epd (= x 1.0d0))
         (x (+ y (* d (/ eps 2.0d0))))
         (epd1 (not (= x y))))
    (when (not (equalp epd epd1))
      (pnl "Something about the next statement is wrong!"))
    (let ((epds (if epd
                    "discarded from most (sub)expressions"
                  (let* ((x (expt 0.5d0 1022))
                         (z (/ (* x (- 2.0d0 (* 2.0d0 eps))) 
                               (- 2.0d0 eps))))
                    (if (< z x)
                        "turned off or do not exist"
                      "turned off")
                    ))))
      (princ "Extra-precise digits are ")
      (princ epds)
      (terpri))
    ))
         
(defun divulge ()
  (progn
    (t1)
    (t2)
    (t3)
    (t4)
    (t5)
    (t6)
    (values)))
