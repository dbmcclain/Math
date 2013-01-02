
(defconstant $gammln-coffs
  #(76.18009173D0
    -86.50532033D0
    24.01409822D0
    -1.231739516D0
    0.120858003D-2
    -0.536382D-5))

(defun gammln (x)
  ;; log Gamma for x > 0
  ;; for x > 1 we use a power series and gain full accuracy
  ;; for 0 < x < 1 the reflection formula is used:  Gamma(1-z) = Pi*z / (Gamma(1+z)*sin(Pi*z))
  (if (< x 1.0d0)
      (let ((z (- 1.0d0 x)))
        (/ (* pi z)
           (* (gammln (+ 1.0d0 z)) (sin (* pi z)))
           ))
    (let* ((xm1 (- x 1.0d0))
           (tmp (+ xm1 5.5d0))
           (tmp (- tmp (* (+ xm1 0.5d0) (log tmp))))
           (ser 1.0d0))
      (loop for j from 0 to 5 do
            (incf xm1 1.0d0)
            (incf ser (/ (aref $gammln-coffs j) xm1)))
      (- (log (* ser 2.50662827465D0)) tmp)
      )))

(defvar *factrl-memo-top*  nil)
(defvar *factrl-memo-cache*
  (let ((cache (make-array 33 :initial-element 0.0d0)))
    (loop for ix from 0
          for val in '(1.0d0 1.0d0 2.0d0 6.0d0 24.0d0)
          do  (setf (aref cache ix)   val
                    *factrl-memo-top* ix))
    cache))

(defun neg-factorial-error (routine)
  (error "Negative factorial in routine ~A" routine))

(defun factrl (n)
  (cond ((minusp n) (neg-factorial-error "FACTRL"))
        ((> n 32)   (exp (gammln (+ 1.0d0 n))))
        (t          (loop
                     while (< *factrl-memo-top* n)
                     do
                     (let ((j *factrl-memo-top*))
                       (incf *factrl-memo-top*)
                       (setf (aref *factrl-memo-cache* *factrl-memo-top*)
                             (* *factrl-memo-top* (aref *factrl-memo-cache* j)))
                       ))
                    (aref *factrl-memo-cache* n))
        ))

(defvar *factln-cache*
  (make-array 101 :initial-element 0.0d0))

(defun factln (n)
  (cond ((minusp n)    (neg-factorial-error "FACTLN"))
        ((<= n 1)      0.0d0)
        ((<= n 100)    (let ((ans (aref *factln-cache* n)))
                         (if (zerop ans)
                             (setf (aref *factln-cache* n) (gammln (+ 1.0d0 n)))
                           ans)))
        (t              (gammln (+ 1.0d0 n)))
        ))

(defun bico (n k)
  ;; binomial coefficient = n! / (k! (n-k)!)
  (values (floor (+ 0.5d0
                    (exp (- (factln n)
                            (factln k)
                            (factln (- n k)))
                         ))
                 )))

(defun beta (z w)
  (exp (- (+ (gammln z) (gammln w))
          (gammln (+ z w)))
       ))


