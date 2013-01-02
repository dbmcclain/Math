
(defconstant $sqrt2  (sqrt 2d0))
(defconstant $sqrtpi (sqrt pi))
(defconstant $root5-double-float-epsilon (* (sqrt 5d0) double-float-epsilon))

(defun lngamma (x)
  
(defun gamma-inc-p (a x)
  (assert (plusp a))
  (assert (not (minusp x)))
  
  (cond ((zerop x)  0d0)

        ((or (< x 20d0)
             (< x (* 0.5d0 a)))
         ;; easy series cases. Robust and quick.
         (gamma-inc-p-series a x))

        ((and (> a 1d6)
              (< (sq (- x a)) a))
         ;; crossover region. Note that Q and P are
         ;; roughly the same order of magnitude here,
         ;; so the subtraction is stable.
         (- 1d0 (gamma-inc-q-asymp-unif a x)))

        ((<= a x)
         ;; Q <~ P in this area, so the
         ;; subtractions are stable.
         (- 1d0
            (if (> a (* 0.2d0 x))
                (gamma-inc-q-cf a x)
              (gamma-inc-q-large-x a x))))

        ((< (sq (- x a)) a)
         ;; This condition is meant to ensure
         ;; that Q is not very close to 1,
         ;; so the subtraction is stable.
         (- 1d0 (gamma-inc-q-cf a x)))

        (t
         (gamma-inc-p-series a x))
        ))

(defun gamma-inc-q-cf (a x)
  (let ((d (gamma-inc-d a x))
        (f (gamma-inc-f-cf a x)))
    (* d (/ a x) f)
    ))

(defun gamma-inc-d (a x)
  ;; the dominant part,
  ;; D(a,x) := x^a e^(-x) / Gamma(a+1)
  (cond ((< a 10d0)
         (exp (- (* a (log x))
                 x
                 (lngamma (+ 1d0 a)))))

        (t
         (let* ((ln-term (if (< x (* 0.5d0 a))
                             (let* ((u    (/ x a))
                                    (ln-u (log u)))
                               (+ 1d0 (- ln-u u)))
                           (let ((mu (/ (- x a) a)))
                             (log-1plusx-mx mu) ;; log(1+mu) - mu
                             )))
                (gstar (gammastar a))
                (term1 (/ (exp (* a ln-term))
                          (sqrt (* 2d0 pi a)))))
           (/ term1 gstar)
           ))
        ))

(defun gamma-inc-p-series (a x)
  (let* ((nmax 5000)
         (d (gamma-inc-d a x)))
    (do* ((n    0   (1+ n))
          (term 1d0 (* term (/ x (+ a n))))
          (sum  1d0 (+ sum term)))
         ((< (abs (/ term sum)) double-float-epsilon)
          (* d sum))
      (assert (< n nmax))
      )))

(defun gamma-inc-q-large-x (a x)
  ;; Q large x asymptotic
  (let ((nmax 5000)
        (d  (gamma-inc-d a x)))
    (do* ((n 0 (1+ n))
          (last 1d0 term)
          (term 1d0 (* term (/ (- a n) x)))
          (sum  1d0))
         ((or (> (abs (/ term last)) 1d0)
              (< (abs (/ term sum)) double-float-epsilon))
          (* d (/ a x) sum))
      (assert (< n nmax))
      (incf sum term))
    ))

(defun gamma-inc-q-asymp-unif (a x)
  ;; Uniform asymptotic for x near a, a and x large.
  ;; See [Temme, p. 285]
  (let* ((rta (sqrt a))
         (eps (/ (- x a) a))
         (ln-term (log-1plusx-mx eps)) ;; log(1+eps) - eps
         (eta (* (signum eps) (sqrt (- (* 2d0 ln-term)))))
         (erfc (erfc (/ (* eta rta) $sqrt2))))
    (multiple-value-bind (c0 c1)
        (cond ((< (abs eps) $root5-double-float-epsilon)
               (values (+ (/ -3d0)
                          (* eps
                             (- (/ 12d0)
                                (* eps (- (/ 23d0 540d0)
                                          (* eps (- (/ 353d0 12960d0)
                                                    (* eps (/ 589d0 30240d0))
                                                    ))
                                          ))
                                )))
                       (- (/ -540d0) (/ eps 288d0))))
              (t
               (let ((rt-term (sqrt (/ (* -2d0 ln-term) (* eps eps))))
                     (lam (/ x a)))
                 (values (/ (- 1d0 (/ 1d0 rt-term)) eps)
                         (/ (- (- (* eta eta eta
                                     (+ (* lam lam)
                                        (* 10d0 lam)
                                        1d0))
                                  (* 12d0 eps eps eps)))
                            (* 12d0 eta eta eta)))
                 )))
      (+ (* 0.5d0 erfc) (* (/ (exp (* -0.5d0 a eta eta))
                              (* $sqrt2 $sqrtpi rta))
                           (+ c0 (/ c1 a))))
      )))

(defun gamma-inc-f-cf (a x)
  ;; Continued fraction which occurs in evaluation
  ;; of Q(a,x) or Gamma(a,x).
  ;;
  ;;              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
  ;;   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
  ;;             1 +   1 +     1 +   1 +      1 +   1 +
  ;;
  ;; Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
  ;;
  ;; Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
  ;; See gamma_inc_Q_CF() below.
  (let* ((nmax  5000)
         (small (expt double-float-epsilon 3))
         (hn    1d0)
         (cn    (/ small))
         (dn    1d0))
    (do  ((n 2 (1+ n))
          (done nil))
        (done hn)
      (assert (< n nmax))
      (let ((an (if (oddp n)
                    (/ (* 0.5d0 (1- n)) x)
                  (/ (- (* 0.5d0 n) a) x))))
        (setf dn (let ((dn (+ 1d0 (* an dn))))
                   (if (< (abs dn) small)
                       small
                     dn)))
        (setf cn (let ((cn (+ 1d0 (/ an cn))))
                   (if (< (abs cn) small)
                       small
                     cn)))
        (setf dn (/ dn))
        (let ((delta (* cn dn)))
          (setf hn (* hn delta))
          (setf done (< (abs (- delta 1d0)) double-float-epsilon)))
        ))
    ))

(defun gamma-inc-q-cf (a x)
  ;; Continued fraction for Q.
  ;;
  ;; Q(a,x) = D(a,x) a/x F(a,x)
  ;;
  ;; Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
  ;;
  ;; Since the Gautschi equivalent series method for CF evaluation may lead
  ;; to singularities, I have replaced it with the modified Lentz algorithm
  ;; given in
  ;;
  ;; I J Thompson and A R Barnett
  ;; Coulomb and Bessel Functions of Complex Arguments and Order
  ;; J Computational Physics 64:490-509 (1986)
  ;;
  ;; In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
  ;; removed.
  ;;
  ;; Identification of terms between the above equation for F(a, x) and
  ;; the first equation in the appendix of Thompson&Barnett is as follows:
  ;;
  ;;    b_0 = 0, b_n = 1 for all n > 0
  ;;
  ;;    a_1 = 1
  ;;    a_n = (n/2-a)/x    for n even
  ;;    a_n = (n-1)/(2x)   for n odd
  ;;
  (let* ((d (gamma-inc-d a x))
         (f (gamma-inc-f-cf a x)))
    (* d (/ a x) f)
    ))

(defun gamma-inc-q-series (a x)
  ;; Useful for small a and x. Handle subtraction analytically.
  (assert nil))


            
