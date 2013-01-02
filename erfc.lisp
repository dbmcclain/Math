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
;; Caller supplies a function fnxy that, when furnished with the index ix and a state,
;; returns the next numerator x[ix] and denominator y[ix], for index ix = 1, 2, ...,
;; along with an updated state.
;;
;; Also required are the initial-state, and the starting values x[0] and y[0].
;;
;; Iteration ceases when two successive iterations produce the same answer.
;; That answer is supplied as the result.
;;
(defun cfrac-iter (fnxy initial-state x0 y0)
  ;; fnxy is a function that should take index ix = 1,2,... and some state,
  ;; and should return a list of x[ix], y[ix], and the new-state.
  (um:with-tail-pure-code
    (labels ((iter (ix state p1 q1 p2 q2)
               (destructuring-bind (x0 y0 new-state) (funcall fnxy ix state)
                 (let* ((p0 (+ (* y0 p1) (* x0 p2)))
                        (q0 (+ (* y0 q1) (* x0 q2)))
                        (v  (/ p0 q0)))
                   (if (= v (/ p1 q1))
                       v
                     (iter (1+ ix) new-state p0 q0 p1 q1))
                   ))))
      (iter 1 initial-state x0 y0 0.0d0 1.0d0)
      )))


;; APriori Continued Fractions for Erf(x) and Erfc(x)
(defun erf_raw (x)
  ;; use when x <= 2.75
  (let* ((xsq (* x x))
         (fnxy (lambda (ix sgn)
                 ;; state is a sign that alternates +1,-1
                 (list (* sgn -2 ix xsq)  ;; x[ix] = numerator for this iter
                       (+ ix ix 1)        ;; y[ix] = denominator for this iter
                       (- sgn)))))        ;; new state is negated incoming sign integer
    (* (/ 2.0d0 (sqrt pi))
       (exp (- xsq))
       (cfrac-iter fnxy 1 x 1.0d0))
    ))

(defun erfc_raw (x)
  ;; use when x > 2.75
  (let ((fnxy (lambda (ix state)
                ;; state unused here
                (list (* 0.5d0 ix)       ;; x[ix] = numerator for this iter
                      x                  ;; y[ix] = denominator for this iter
                      state))))          ;; just return the unused state
    (* (/ (sqrt pi))
       (exp (- (* x x)))
       (cfrac-iter fnxy nil 1.0d0 x))
    ))

;; -------------------------------------------------------------
;; User callable entry points
;;
;; These entry points determine which of the two raw definitions to call,
;; based on the magnitude of the argument x. When abs(x) = 2.75 both raw
;; definitions require about the same number of iterations for convergence.
;;
(defun erfc (x)
  ;; 2/Sqrt(Pi)*Integral(Exp(-t^2), {t, x, inf}) = 1 - erf(x)
  (let* ((z   (abs (float x 1.0d0)))
        (ans  (if (> z 1.0d0)
                  ;; using a boundary of 1.0 instead of 2.75 gives about
                  ;; equal hair on each side of the boundary
                  (erfc_raw z)
                (- 1.0d0 (erf_raw z)))))
    (if (minusp x)
        (- 2.0d0 ans)
      ans)))

(defun erf (x)
  ;; 2/Sqrt(Pi)*Integral(Exp(-t^2), {t, 0, x}) = 1 - erfc(x)
  (let* ((z   (abs (float x 1.0d0)))
         (ans (if (> z 1.0d0)
                  ;; using a boundary of 1.0 instead of 2.75 gives about
                  ;; equal hair on each side of the boundary
                  (- 1.0d0 (erfc_raw z))
                (erf_raw z))))
    (if (minusp x)
        (- ans)
      ans)))

#| ;check it out...
(let ((domain '(-3.0d0 3.0d0)))
  (plt:fplot 'erfc1 domain #'erfc
             :title "Erfc(x)")
  
  (plt:fplot 'erfc2 domain #'erf
            :title "Erf(x)")
  
  (plt:fplot 'erfc3 domain (lambda (x)
                     (- (stocks::erfc x)
                        (erfc x)))
            :title "Erfc Approximation Error"))

|#