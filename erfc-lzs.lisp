;; Erfc.lisp -- Apriori canonical routines for Erf(x) and Erfc(x)
;; These are the definitional infinite continued fraction representations for the
;; integrals. These routines may be used to ascertain the errors in faster approximations.
;;
;; DM/RAL 12/06
;; ----------------------------------------------------------------------

;; ----------------------------------------------------------------------
;; cfrac-iter -- an evaluator for continued fraction approximations,
;; expressed in the form
;;
;;  f(x) = x0
;;         -------
;;         y0 + x1
;;              -------
;;              y1 + x2
;;                   -------
;;                   y2 + x3 .....
;;
;; could also be written as:
;;
;; f(x) = x0/(y0 + x1/(y1 + x2/(y2 + ...
;;
;; Caller supplies a function fnxy that, when furnished with the index ix,
;; returns the next numerator x[ix] and denominator y[ix], for index ix = 1, 2, ...
;;
;; Also required are the starting values x[0] and y[0].
;;
;; We use accelerated iteration with Euler's method.
;; Iteration ceases when two successive iterations produce the same answer.
;; That answer is supplied as the result.
;;
;; This is made easier by the use of lazy-streams.
;;

(defun cfrac-term (fnxy)
  (lambda (args)
    (destructuring-bind (ix v p1 q1 p2 q2) args
      (declare (ignore v))
      (destructuring-bind (x0 y0) (funcall fnxy ix)
        (let ((p0 (+ (* y0 p1) (* x0 p2)))
              (q0 (+ (* y0 q1) (* x0 q2))))
          (list (1+ ix) (/ p0 q0) p0 q0 p1 q1)
          )))))

(defun repeat (f a)
  (lzs:stream-cons a (repeat f (funcall f a))))

(defun erf-stream (x)
  ;; cfrac is: x|1-2*x^2|3+4*x^2|5-6*x^2|7+ ...
  ;; use when x <= 1.7
  (let* ((xsq (* x x))
         (fnxy (lambda (ix)
                 (let ((sgn (- 1 (* 2 (logand ix 1)))))
                   (list (* sgn 2 ix xsq) (+ 1 ix ix))
                   ))))
    (lzs:stream-map #'second (repeat (cfrac-term fnxy) (list 1 x x 1 0 1)))
    ))

(defun erfc-stream (x)
  ;; cfrac is: 1|x+(1/2)|x+(2/2)|x+(3/2)|x+ ...
  ;; use when x > 1.7
  (let ((fnxy (lambda (ix)
                (list (/ ix 2) x))))
    (lzs:stream-map #'second (repeat (cfrac-term fnxy) (list 1 (/ x) 1 x 0 1)))
    ))

(defun erf-raw (x)
  ;; use when x <= 1.7
  (* (/ 2.0d0 (sqrt pi) (exp (* x x)))
     (lzs:within #'=
                 (lzs:accelerate-series #'lzs:euler
                                        (erf-stream x))
                 )))

(defun erfc-raw (x)
  ;; use when x > 1.7
  (* (/ 1.0d0 (sqrt pi) (exp (* x x)))
     (lzs:within #'=
                 (lzs:accelerate-series #'lzs:euler
                                        (erfc-stream x))
                 )))

;; -------------------------------------------------------------
;; User callable entry points
;;
;; These entry points determine which of the two raw definitions to call,
;; based on the magnitude of the argument x. When abs(x) = 1.7 both raw
;; definitions require about the same number of accelerated iterations for convergence.
;;
(defun erfc (x)
  ;; 2/Sqrt(Pi)*Integral(Exp(-t^2), {t, x, inf}) = 1 - erf(x)
  (let* ((z   (abs (float x 1.0d0)))
        (ans  (if (> z 1.7d0)
                  (erfc-raw z)
                (- 1.0d0 (erf-raw z)))))
    (if (minusp x)
        (- 2.0d0 ans)
      ans)))

(defun erf (x)
  ;; 2/Sqrt(Pi)*Integral(Exp(-t^2), {t, 0, x}) = 1 - erfc(x)
  (let* ((z   (abs (float x 1.0d0)))
         (ans (if (> z 1.7d0)
                  (- 1.0d0 (erfc-raw z))
                (erf-raw z))))
    (if (minusp x)
        (- ans)
      ans)))

#| ;check it out...
(let ((domain '(-3.0d0 3.0d0)))
  (plt:fplot 1 domain #'erfc
             :clear t
             :title "Erfc(x)")
  
  (plt:fplot 2 domain #'erf
             :clear t
             :title "Erf(x)")
  
  (plt:fplot 3 domain (lambda (x)
                        (- (stocks::erfc x)
                           (erfc x)))
             :clear t
             :title "Erfc Approximation Error"))

|#