;; primes.lisp -- generation and testing of prime numbers
;; --------------------------------------------------------------------------------------
;; DAM Files -- direct access disk-based hash tables.
;;
;; Copyright (C) 2008 by SpectroDynamics, LLC. All rights reserved.
;;
;; DM/SD  08/08
;; --------------------------------------------------------------------------------------

(defpackage :primes
  (:use #:common-lisp)
  (:export
   #:divides?
   #:expt-mod
   #:random-between
   #:make-prime
   #:is-prime?
   #:extended-gcd
   #:compute-modulo-inverse
   #:provably-prime?
   #:factors-of
   #:generate-strong-prime
   #:generate-rsa-base
   #:expt-mod
   #:inv-mod
   #:mult-mod
   ))

;; -------------------------------------------
(in-package #:primes)
;; -------------------------------------------

(defun mt-random (ix)
  (declare (integer ix))
  ;; return a random integer 0 <= val < ix
  #+:LISPWORKS (lw:mt-random ix)
  #+:ALLEGRO   (random ix)) ;; Allegro uses Mersenne Twister already

(defun random-between (lower upper)
  ;; generate random integer in [lower, upper)
  (+ lower (mt-random (- upper lower))))

(defun divides? (n d)
  (declare (integer n d))
  (and (<= d n)
       (zerop (mod n d))))

(defun expt-mod (base exponent modulus)
  ;; base^exponent mod modulus, for any modulus
  #F
  (declare (integer base exponent modulus))
  (let ((n (integer-length exponent)))
    (declare (integer n))
    (do ((b  base (mult-mod b b modulus))
         (p  1)
         (ix 0    (1+ ix)))
        ((>= ix n) p)
      (declare (integer b p ix))
      (when (logbitp ix exponent)
        (setf p (mult-mod p b modulus)))) ))
#|
(defun expt-mod (a e m)
  ;; (a^e) mod m
  (declare (integer a e m))
  (do ((ix 1 (ash ix 1))
       (prod 1)
       (mult a (mod (* mult mult) m)))
      ((> ix e) prod)
    (declare (integer ix prod mult))
    (unless (zerop (logand ix e))
      (setf prod (mod (* mult prod) m)))
    ))
|#

(defun mult-mod (a b m)
  (mod (* a b) m))

(defun inv-mod (a n)
  ;; 1/a mod n
  ;; assumes n is prime
  ;; solve Bezouts identity: a*x + b*y = gcd(a,b)
  ;; for (x, y) => x is the modular multiplicative inverse of a modulo b,
  ;; when a coprime to b.
  #F
  (declare (integer a n))
  (unless (< a n)
    (error "INV-MOD A not less than Modulus"))
  (let ((i  n)
        (j  a)
        (y1 1)
        (y2 0))
    (declare (integer i j y2 y1))
    (loop for (quotient remainder) of-type (integer integer)
                    = (multiple-value-list (truncate i j))
          do
          (when (zerop remainder)
            (unless (= j 1)
              (error "INV-MOD: no inverse exists"))
            (return (mod y1 n)))
          (shiftf i j remainder)
          (shiftf y2 y1 (- y2 (* y1 quotient)) )) ))

;; ----------------------------------------------------------

#|
(defvar *primes*
  ;; first 168 primes, good for testing up to 1,000,000.
  (do ((s   '(19 17 13 11 7 5 3 2))
       (ix   23 (+ ix 2)))
      ((>= ix 1000) (nreverse s))
    (unless (member ix s
                    :test #'divides?)
      (push ix s))
    ))
|#

(defvar *primes*
  ;; first 168 primes, good for testing up to 1,000,000.
  (do ((s   '(3 2)) ;; all other primes are (6*k +/- 1)
       (ix   6 (+ 6 ix)))
      ((>= ix 1000) (nreverse s))
    (let ((p (1- ix)))
      (unless (member p s
                      :test #'divides?)
      (push p s)))
    (let ((q (1+ ix)))
      (unless (member q s
                      :test #'divides?)
        (push q s))) ))

;; ----------------------------------------------------------

(defun make-prime (k &optional (mr-iters 50))
  ;; return a prime p, k < p < 2*k-2
  ;; By Bertrand's Postulate, there is always at least one prime, p,
  ;; where k < p < 2*k-2, for k > 3
  (declare (integer k))
  (cond ((< k 997)
         (let ((ps (reduce (lambda (ps p)
                             (if (< k p (* 2 (1- k)))
                                 (cons p ps)
                               ps))
                           *primes*
                           :initial-value nil)))
           (values (elt ps (mt-random (length ps))) 1)))

        (t
         (let ((lower (1+ k))
               (upper (* 2 (1- k)))
               (niter 0))
           
           ;; Since fundamental theorem states that pi(n) ~ n/ln(n),
           ;; or density of primes is inversely proportional to number of decimal digits,
           ;; on average about ln(k)/2 trials will be needed
           #| |#
           (let ((n (1- (* 6 (ceiling (random-between lower upper) 6)))))
             (values (or (um:nlet-tail iter> ((n n))
                           (incf niter)
                           (when (< n upper)
                             (cond ((is-prime? n mr-iters) n)
                                   ((is-prime? (+ n 2) mr-iters)
                                    (incf niter)
                                    (+ n 2))
                                   (t
                                    (incf niter)
                                    (iter> (+ n 6))) )))
                         (um:nlet-tail iter< ((n (- n 6)))
                           (incf niter)
                           (when (> n lower)
                             (cond ((is-prime? n mr-iters) n)
                                   ((is-prime? (+ n 2) mr-iters)
                                    (incf niter)
                                    (+ n 2))
                                   (t
                                    (incf niter)
                                    (iter< (- n 6))) ))))
                     niter))
           #| |#
           #|
           (let ((n (logior 1 (random-between lower upper))))
             (values (or (um:nlet-tail iter> ((n n))
                           (incf niter)
                           (when (< n upper)
                             (cond ((is-prime? n) n)
                                   (t (iter> (+ n 2))) )))
                         (um:nlet-tail iter< ((n (- n 2)))
                           (incf niter)
                           (when (> n lower)
                             (cond ((is-prime? n) n)
                                   (t (iter< (- n 2))) ))))
                     niter))
           |#
           ))
        ))

#|
(defun collect-stats (n k)
  (let ((trials (loop repeat n collect (second (multiple-value-list (make-prime k))))))
    (plt:histogram 'plt2 trials
                   :clear t
                   :thick 2
                   :title  "Number of Probes to Find a 100-digit Prime"
                   :ytitle "PDF Density"
                   :xtitle "Number of PRNG Probes")))
|#

;; ----------------------------------------------------------
#|
(defun is-prime? (p)
  (let ((p-1 (1- p)))
    (labels ((witness? (witness)
               (/= 1 (expt-mod witness p-1 p))))
      
      (cond ((<= p 1000000)
             (or (member p *primes*)
                 (not (member p *primes*
                              :test #'divides?))
                 ))

            ((evenp p) nil) ;; no higher primes are even
          
            ((zerop (mod p 5)) nil) ;; no higher primes end in 5
        
            ((some #'witness? *primes*) nil) ;; lower primes are highly likely witnesses
            
            (t
             ;; use probabilistic Fermat's test for certainty < 1/2^100
             ;; if p prime then a^(p-1) mod p = 1, for 1 <= a < p
             ;;
             ;; NOTE: this fails if p is a Charmichael number = a composite
             ;; number for which a^(p-1) mod p = 1 whenever gcd(a,p) = 1
             (um:nlet-tail iter ((ntests 100)
                                 (seen   *primes*))
               (if (zerop ntests)
                   t ;; maybe...
                 (let ((a (+ 4 (mt-random (- p 5)))))
                   (if (member a seen)
                       (iter ntests seen)
                     (if (witness? a)
                         nil
                       (iter (1- ntests) (cons a seen))) )))))
            ))))
                              
;; --------------------------------------------------------------

(defun jacobi-symbol (m n &optional (prod 1))
  ;; requires n odd
  (assert (oddp n))
  (if (> m n)
      (jacobi-symbol (mod m n) n prod)
    
    (cond
     ((= 1 n) prod)
     
     ((= 1 m)
      (if (logbitp 1 n)
          (- prod)
        prod))
     
     ((= 2 m)
      (if (member (mod n 8) '(1 7))
          prod
        (- prod)))
     
     ((evenp m)
      (jacobi-symbol (ash m -1) n (jacobi-symbol 2 n prod)))
     
     ((/= 1 (gcd m n)) 0)
     
     (t 
      (jacobi-symbol n m
                     (if (and (= 3 (mod n 4))
                              (= 3 (mod m 4)))
                         (- prod)
                       prod)))
     )))

;; --------------------------------------------------------------------

(defun solovay-strassen-prime? (p &optional (k 100))
  ;; false positive < 1/2^k
  ;; require p odd > 2
  (assert (oddp p))
  (assert (> p 2))
  (um:nlet-tail iter ((niter k))
    (if (zerop niter)
        t
      (let ((a (+ 2 (mt-random (- p 2)))))
        (when (= 1 (gcd a p))
          (let ((jac  (jacobi-symbol a p))
                (pwr  (expt-mod a (ash p -1) p)))
            (when (= jac pwr)
              (iter (1- niter))) ))) )))
|#

#|
(defun random-in-range (min max)
  "Return a random integer in the interval [min .. max]"
  (declare (integer min max))
  (+ min (mt-random (1+ (- max min)))))
|#

#|
(defun make-unique-random-in-range-fn (min max)
  (declare (integer min max))
  (let ((seen  nil)
        (limit (1+ (- max min))))
    (declare (integer limit))
    (lambda ()
      (um:nlet-tail iter ()
        (let ((n (mt-random limit)))
          (declare (integer n))
          (if (member n seen :test #'=)
              (iter)
            (progn
              (push n seen)
              (+ min n))) ))) ))
|#

(defun factor-out (number divisor)
  "Return two values R and E such that NUMBER = DIVISOR^E * R,
  and R is not divisible by DIVISOR."
  (declare (integer number divisor))
  (do ((e 0 (1+ e))
       (r number (/ r divisor)))
      ((/= (mod r divisor) 0) (values r e))
    (declare (integer e r)) ))

(defun perfect-square? (c)
  (let* ((n  (integer-length c))
         (m  (ceiling n 2)))
    (do ((x  (random-between (ash 1 (1- m)) (ash 1 m))
             (truncate (+ c (* x x)) (* 2 x))))
        ((< (* x x) (+ c (ash 1 m)))
         (values (= c (* x x)) x))
      )))

(defun jacobi-symbol (a n)
  (let ((a (mod a n)))
    (cond
     ((or (= 1 a)
          (= 1 n))
      1)
     
     ((zerop a) 0)
     
     (t (multiple-value-bind (a1 e) (factor-out a 2)
          (let* ((rem (mod n 8))
                 (s   (cond
                       ((evenp e) 1)
                       ((or (= 1 rem)
                            (= 7 rem)) 1)
                       ((or (= 3 rem)
                            (= 5 rem)) -1))))
            (when (and (= 3 (mod n 4))
                       (= 3 (mod a1 4)))
              (setf s (- s)))
            (* s (jacobi-symbol (mod n a1) a1)) )))
     )))
                 
(defun probabilistic-lucas-test (c)
  (if (perfect-square? c)
      nil
    (let* ((d (loop for sign = 1 then (- sign)
                    for d = 5 then (* sign (+ 2 (abs d)))
                    for jac = (jacobi-symbol d c)
                    do
                    (cond ((zerop jac) (return-from probabilistic-lucas-test nil))
                          ((= -1 jac)  (return d)))))
           (k  (1+ c)))

      (labels ((half-mod (x)
                 (mod (ash (if (oddp x)
                               (+ x c)
                             x)
                           -1)
                      c)))
        
        (um:nlet-tail iter ((u  1)
                            (v  1)
                            (ix (- (integer-length k) 2)))
          (if (minusp ix)
              (zerop u)
            (let ((utmp (mult-mod u v c))
                  (vtmp (half-mod (+ (* v v) (* d u u)))))
              (if (logbitp ix k)
                  (iter (half-mod (+ utmp vtmp))
                        (half-mod (+ vtmp (* d utmp)))
                        (1- ix))
                (iter utmp vtmp (1- ix))) ))) ))))
                
      
#|
(defun strong-miller-rabin-liar? (d s p p-1 a)
  ;; return true if a not not witness a composite value for p
  ;; p-1 = d * 2^s
  (let ((x (expt-mod a d p)))
    (or (= x 1)
        (loop repeat s
              for y = x then (mod (* y y) p)
              thereis (= y p-1))) ))
|#

(defun make-strong-miller-rabin-liar-test (p)
  (declare (integer p))
  (let ((p-1 (1- p)))
    (declare (integer p-1))
    (multiple-value-bind (d s) (factor-out p-1 2)
      ;; p-1 = d * 2^s, d odd
      (declare (integer d s))
      
      (lambda (a)
        ;; return true if a not not witness a composite value for p
        (declare (integer a))
        (let ((x (expt-mod a d p)))
          (declare (integer x))
          (or (= x 1)
              (loop repeat s
                    for y of-type integer = x then (the integer (mod (the integer (* y y)) p))
                    thereis (= y p-1))) )))))

(defun probabilistic-miller-rabin-prime? (p &optional (k (min 50 (- p 2))))
  ;; false positive < 1/4^k
  ;; return NIL if p found to be composite
  ;; assumes p odd > 2
  (declare (integer p k))
  (let ((p-1 (1- p)))
    (declare (integer p-1))
    (assert (> p-1 k))
    (let ((strong-liar? (make-strong-miller-rabin-liar-test p)))
      (loop repeat k
            for witness of-type integer = (random-between 2 p)
            always (funcall strong-liar? witness)) )))

(defun deterministic-miller-rabin-prime? (p witnesses)
  (declare (integer p))
  ;; return NIL if p found to be composite on the basis of the listed witnesses
  ;; assumes p odd > 2
  (let ((strong-liar? (make-strong-miller-rabin-liar-test p)))
    (loop for witness of-type integer in witnesses
          always (funcall strong-liar? witness)) ))

(defun miller-rabin-prime? (p &optional (mr-iters 50))
  ;; assumes p odd > 3
  (declare (integer p))
  (cond
   ((< p         1373653) (deterministic-miller-rabin-prime? p '(2 3)))
   ((< p         9080191) (deterministic-miller-rabin-prime? p '(31 73)))
   ((< p      4759123141) (deterministic-miller-rabin-prime? p '(2 7 61)))
   ((< p   2152302898747) (deterministic-miller-rabin-prime? p '(2 3 5 7 11)))
   ((< p   3474749660383) (deterministic-miller-rabin-prime? p '(2 3 5 7 11 13)))
   ((< p 341550071728321) (deterministic-miller-rabin-prime? p '(2 3 5 7 11 13 17)))
   (t                     (and (probabilistic-miller-rabin-prime? p mr-iters)
                               (probabilistic-lucas-test p)))
   ))


(defun is-prime? (p &optional (mr-iters 50))
  (declare (integer p))
  (cond
   ((= p 1)   nil)
   ((< p 4)   t)
   ((evenp p) nil)
   (t         (miller-rabin-prime? p mr-iters))
   ))

;; ------------------------------------------------------------

(defun extended-gcd (a b)
  ;; solve Bezout's identity: a*x + b*y = gcd(a,b)
  ;; for (x, y) => x is the modular multiplicative inverse of a modulo b,
  ;; when a coprime to b.
  (if (zerop (mod a b))
      (values 0 1)
    (multiple-value-bind (x y) (extended-gcd b (mod a b))
      (values y
              (- x (* y (truncate a b)))) )))

(defun compute-modulo-inverse (a m)
  (assert (= 1 (gcd a m))) ;; ensure A coprime with M
  (mod (extended-gcd a m) m))

;; ---------------------------------------------------------------
#|
(defun tst (n)
  (let ((u (expt 65537 10)))
    (loop repeat n
          for ix from 1
          do
          (let ((y (1+ (* 2 ix u))))
            (when (is-prime? y)
              (print y)
              (force-output))))))
|#
;; ---------------------------------------------------------------
;; (defvar *niter* 0)

(defun test-primitive-root (n factors x)
  ;; check to see if x is a primitive root of prime n
  ;; where factors are the prime factors of phi(n) = (n-1)
  (and
   ;; (= 1 (gcd x n))                  ;; - coprime test S.B. true for all X - REDUNDANT
   (= 1 (expt-mod x (1- n) n))      ;; - Fermat's test S.B. true for all X
   (every (lambda (factor) ;; - Lucas' test S.B. true for some X
            (/= 1 (expt-mod x (truncate (1- n) (first factor)) n)))
          factors)
   x)) ;; return the witness as true
  
(defun lucas-lehmer-test (n factors &optional (limit 20))
  ;; Lucas-Lehmer test for N = prod(factors)+1 as prime,
  ;; for FACTORS a list of prime factors of (N-1).
  ;; If passes test, then N is absolutely prime and we return a witness number.
  ;; Otherwise it is probably composite and we return NIL.
  ;;
  ;; This test can give a false positive on composite numbers if
  ;; the list of factors of (N-1) include composite numbers.
  ;; So be certain that factors are all prime numbers.

  ;; (incf *niter*)
  (or ;; (test-primitive-root n factors 2)
      ;; (test-primitive-root n factors 5)
      (um:nlet-tail iter ((limit limit))
        (and (plusp limit)                    ;; - true if not exhausted the search
             (let ((x (random-between 2 (1- n))))
               ;; no need to look over whole range, since the upper half mirrors the lower half
               (or (test-primitive-root n factors x) ;; - Lucas' test S.B. true for some X
                   (iter (1- limit)))) ))))

(defun provably-prime? (n &optional absolute-proof)
  (lucas-lehmer-test n (factors-of (1- n) absolute-proof)))

#|
(defun collect-stats (n k)
  (let ((trials (loop repeat n collect
                      (let ((*niter* 0))
                        (provably-prime? (make-prime k))
                        *niter*)) ))
    (plt:histogram 'plt3 trials
                   :clear t
                   :thick 2
                   :cum    nil
                   :title  "Number of Probes in Lucas-Lehmer Test"
                   :ytitle "PDF Density"
                   :xtitle "Number of PRNG Probes")
    (values (float (vm:mean trials))
            (float (vm:stdev trials))) ))


(defun collect-stats (n k)
  (labels ((hit-test (&optional (niter 1))
               (if (divides? (mt-random k) 3)
                   niter
                 (hit-test (1+ niter)))))
    (let ((trials (loop repeat n collect
                        (hit-test))))
      (plt:histogram 'plt3 trials
                     :clear t
                     :thick 2
                     :cum    nil
                     :title  "Number of Probes in 1/3 Hit-Test"
                     :ytitle "PDF Density"
                     :xtitle "Number of PRNG Probes")
      (values (float (vm:mean trials))
              (float (vm:stdev trials))) )))
|#

#|
(defun small-composite? (n)
  ;; test for absolute prime N, for N < 1,000,000
  ;; by using direct division tests for primes in *PRIMES*
  (assert (< n 1000000))
  (let ((limit (isqrt n)))
    (some (lambda (x)
            (and (< x limit)
                 (divides? n x)))
          *primes*)))
|#

(defun factors-of (n &optional complete-factorization-p)
  ;; good for some N -- return a list of conses: ((factor1 . exponent1) (factor2 . exponent2) ...)
  ;;
  ;; It is possible that after partial factoring a residue (> 1) will remain.
  ;; 
  ;; If the residue IS NOT a probable prime, then it is definitely composite, and its factors
  ;; are larger than we can handle from *PRIMES*. So give up with an error report.
  ;; 
  ;; If the residue IS a probable prime then:
  ;;   - If the user asked for a complete factorization, attempt to recursively
  ;;       prove that the residue is prime, and then append that residue to the list of factors.
  ;;       If we are unable to prove that the residue is prime then it is too large for this routine.
  ;;
  ;;   - If the user did not ask for a complete factorization then just accept it as a
  ;;     probable prime factor, and append it to the list of factors.
  ;;
  (um:nlet-tail iter ((n  n)
                      (ps *primes*)
                      (factors nil))
    (if ps
        (let ((p (first ps)))
          (if (< n p)
              factors
            (if (divides? n p)
                (multiple-value-bind (d r) (factor-out n p)
                  (iter d (rest ps) (cons (cons p r) factors)))
              (iter n (rest ps) factors))))

      (labels ((residue-too-large ()
                 (error "Residue too large in FACTORS-OF: ~A ~A" n factors)))
    
        (cond
         ((= n 1) factors)
         
         ((is-prime? n)
          (if complete-factorization-p
              ;; user asked for complete factorization. Try if possible.
              (if (provably-prime? n t)
                  (cons (cons n 1) factors)
                ;; else - couldn't prove it, the residue is probably composite
                ;; and it must be too large for this routine
                (residue-too-large))
            ;; else
            (cons (cons n 1) factors)))

         (t (let ((r (isqrt n)))
              (if (and (= n (* r r))
                       (is-prime? r))
                  (cons (cons r 2) factors)
                (residue-too-large)))) ;; non-prime leftover -- routine can't handle it
         )))))

(defun totient (n)
  ;; Euler totient function: phi(n)
  ;; For primitive roots modulo n, there are phi(phi(n)) of them
  ;; For n prime, P, phi(P) = (P-1), and the primitive roots appear
  ;; more or less uniformly distributed over the interval 1 < x < P.
  (reduce (lambda (tot factor)
            (destructuring-bind (p . e) factor
              (* tot (1- p) (expt p (1- e)))))
          (factors-of n)
          :initial-value 1))

#|
(plt:plot 'plt
          *primes*
          (mapcar (lambda (p)
                    (float (/ (totient (1- p)) (1- p))))
                  *primes*)
          :clear t)
|#

;; ---------------------------------------------------------------
#|
(defun generate-strong-prime (ndigits &optional (nlevels 2))
  ;; a strong prime N has a large prime factor in N-1
  ;; NLEVELS indicates how many component (P-1) should be strong
  ;;  NLEVELS: T   = all levels should be strong and proven prime
  ;;           NIL = none should be strong -> probabilistic prime
  ;;           <integer> = nbr of strong component levels - proven or probabilistic result
  ;;                 If the printout shows: Make Base Prime: Absolute
  ;;                    then we will have a proven prime result.
  ;;                 Otherwise, we will have a probabilistic prime result.
  (labels ((find-incremental-prime (u)
             ;; start with large prime U, and compute 2*k*U+1 until prime
             (let ((incr (* 2 u u)))
               (loop for p from (1+ incr) by incr
                     for k from 1
                     until (and (is-prime? p)
                                (provably-prime? p)
                                (< (1- ndigits) (* 0.301 (integer-length p))))
                     ;; return a provably prime number, assuming u were prime
                     finally (return p))))
           
           (base-prime (ndigits)
             (format t "~&Make Base Prime: ")
             (loop for u = (make-prime (random-between (expt 10 ndigits) (* 5 (expt 10 ndigits))))
                   until (handler-case
                             (when (provably-prime? u t)
                               (format t "Absolute") ;; proven prime
                               t) ;; if no error, but not proven prime, then try again
                           (error (cx)
                             (declare (ignore cx))
                             ;; could not completely factor
                             (when (or (null nlevels)
                                       (integerp nlevels))
                               (format t "Probabilistic") ;; accept on probabilistic faith
                               t))) ;; otherwise, return nil to generate another trial prime
                   finally (return u))) )

    ;; get starting large prime
    (let ((u (if (or (< ndigits 32)
                     (and nlevels
                          (integerp nlevels)
                          (< nlevels 2))
                     (null nlevels)) ;; Null levels - accept base prime as U
                 (base-prime (truncate ndigits 2))
               ;; else, recursively generate strong prime for U
               (generate-strong-prime (truncate ndigits 2) (or (and (integerp nlevels)
                                                                    (1- nlevels))
                                                               nlevels))) ))
      (format t "~&Using: ~A" u)
      (find-incremental-prime u)) ))
|#

(defun generate-strong-prime (nbits &optional (nlevels 2))
  ;; a strong prime N has a large prime factor in N-1
  ;; NLEVELS indicates how many component (P-1) should be strong
  ;;  NLEVELS: T   = all levels should be strong and proven prime
  ;;           NIL = none should be strong -> probabilistic prime
  ;;           <integer> = nbr of strong component levels - proven or probabilistic result
  ;;                 If the printout shows: Make Base Prime: Absolute
  ;;                    then we will have a proven prime result.
  ;;                 Otherwise, we will have a probabilistic prime result.
  (labels ((find-incremental-prime (u)
             ;; start with large prime U, and compute 2*k*U+1 until prime
             (let ((incr (* 2 u u)))
               (loop for p from (1+ incr) by incr
                     for k from 1
                     until (and (is-prime? p)
                                (provably-prime? p)
                                (< (1- nbits) (integer-length p)))
                     ;; return a provably prime number, assuming u were prime
                     finally (return p))))
           
           (base-prime (nbits)
             (format t "~&Make Base Prime: ")
             (loop for u = (make-prime (let ((b (ash 1 nbits)))
                                         (random-between b (+ b b))))
                   until (handler-case
                             (when (provably-prime? u t)
                               (format t "Absolute") ;; proven prime
                               t) ;; if no error, but not proven prime, then try again
                           (error (cx)
                             (declare (ignore cx))
                             ;; could not completely factor
                             (when (or (null nlevels)
                                       (integerp nlevels))
                               (format t "Probabilistic") ;; accept on probabilistic faith
                               t))) ;; otherwise, return nil to generate another trial prime
                   finally (return u))) )

    ;; get starting large prime
    (let ((u (if (or (< nbits 107)
                     (and nlevels
                          (integerp nlevels)
                          (< nlevels 2))
                     (null nlevels)) ;; Null levels - accept base prime as U
                 (base-prime (- (truncate nbits 2) 3))
               ;; else, recursively generate strong prime for U
               (generate-strong-prime (- (truncate nbits 2) 3) (or (and (integerp nlevels)
                                                                        (1- nlevels))
                                                                   nlevels))) ))
      (format t "~&Using: (~A) ~A" (integer-length u) u)
      (find-incremental-prime u)) ))

;; ----------------------------------------------------------------------------------           

(defun compute-diffie-hellman (nbits)
  (let* ((p       (generate-strong-prime nbits t))
         (factors (factors-of (1- p)))
         (g       (lucas-lehmer-test p factors))) ;; returns a primitive root if successful
    (if g ;; if lucas-lehmer failed then try again
        (list :DH-P p :DH-G g)
      (compute-diffie-hellman nbits))))

;; ----------------------------------------------------------------------------------           

(defun compute-rsa-d&e (p q)
  (let ((phi    (* (1- p) (1- q)))
        (thresh (integer-length (* p q)))) ;; threshold for e greater than this for good mixing
    (loop for d = (make-prime (1+ (max p q)))
          for e = (inv-mod d phi)
          do
          (when (> e thresh)
            (return (values d e)))
          (when (> d thresh)
            (return (values e d)) ))))

(defun generate-rsa-base (nbits)
  (let* ((nd1 (- (floor nbits 2) 8))
         (p   (generate-strong-prime nd1))
         (nd2 (+ (ceiling nbits 2) 8))
         (q   (generate-strong-prime nd2)))
    (multiple-value-bind (d e)
        (compute-rsa-d&e p q)
      ;; perform trial encryptions/decryptions to verify goodness
      (loop repeat 100 do
            (let ((msg (random-between 2 (* p q))))
              (assert (= msg (expt-mod (expt-mod msg e (* p q)) d (* p q))))))
      (list :N-Key (* p q) ;; Public N-Key
            :E-Key e ;; Public E-Key
            :D-Key d ;; Private D-Key
            :Primes  (list p q) ;; Private prime factors of N-Key
            :Phi-GCD (gcd (1- p) (1- q)) ;; verify that GCD is small
            ))))

(defvar *daves-N-Key* ;; public
  2264622332100115635046544861457968154487858357629655981804341141586114228298294029783844324050989940210772053439566053045896814644059794725879563686363912167075068867395096157661687366233697816122565111)

(defvar *daves-E-Key* ;; public
  246217118187169985272570746886853299629231075233304175666349360972368752967752016255671063892504959532124058155571361367065153880715972807546628534121821786783732922864646413843902219977232097825366741)

;; ------------------------------------

(defvar *daves-D-Key* ;; private
  3336134262755358265069226579401091710115226108761764767081277708392567434516013871621774832824438208861)

(defvar *daves-P-prime* ;; private
  1213152860237679366017990949504570398999532153649335816658883894008526264696809742265586934057363591)

(defvar *daves-Q-prime* ;; private
  1866724636544510727941144178801761797219686359320662768042557145048836731901008435492740832105185576721)

;; --------------------------------------------------------------------------------
#|
(um:nlet-tail iter ((ps *primes*)
                    (coll nil))
  (if (rest ps)
      (let ((p (first ps)))
        (do ((qs (rest ps) (cdr qs)))
            ((null qs))
          (let* ((q (first qs))
                 (n (1+ (* p q))))
            (if (is-prime? n)
                (push n coll)) )))
        (iter (rest ps) coll))
    coll))
|#

;; ----------------------------------------------------------------------------------

(defun generate-safe-prime (nbits)
  ;; a safe prime is of the form (2*U+1) where U is a large prime.
  ;; NLEVELS indicates how many component (P-1) should be strong
  ;;  NLEVELS: T   = all levels should be strong and proven prime
  ;;           NIL = none should be strong -> probabilistic prime
  ;;           <integer> = nbr of strong component levels - proven or probabilistic result
  ;;                 If the printout shows: Make Base Prime: Absolute
  ;;                    then we will have a proven prime result.
  ;;                 Otherwise, we will have a probabilistic prime result.
  (loop for u = (make-prime (ash 1 (- nbits 2)))
        for p = (+ 1 u u)
        for niter from 1
        until (is-prime? p)
        finally (return (list niter p))))

;; -----------------------------------------------------------------------------------

(defun show-group (n)
  (loop for ix from 2 below (1- n) ;; to (truncate (1- n) 2)
        collect
        (loop for k from 1 below n collect (expt-mod ix k n))))


(defun show-gens (n)
  (let ((factors (factors-of (1- n))))
    (loop for ix from 2 below (ceiling n 2)
          for g = (test-primitive-root n factors ix)
          when g collect g)))
        
;; -----------------------------------------------------------------------------------

(defun convert-vector-to-integer (v)
  (let* ((nb (length v))
         (n  0))
    (loop for iy from (* 8 (1- nb)) downto 0 by 8
          for ix from 0
          do
          (setf (ldb (byte 8 iy) n) (aref v ix)))
    n))

(defun convert-integer-to-vector (x)
  (declare (integer x))
  (let* ((nb (ceiling (integer-length x) 8))
         (v  (make-array nb :element-type '(unsigned-byte 8))))
    (declare (fixnum nb))
    (loop for ix of-type fixnum from (* 8 (1- nb)) downto 0 by 8
          for iy of-type fixnum from 0
          do
          (setf (aref v iy) (ldb (byte 8 ix) x)))
    v))

;; -----------------------------------------------

(defmethod integer-of ((x integer))
  x)

(defmethod integer-of ((v vector))
  (convert-vector-to-integer v))

(defmethod integer-of ((u uuid:uuid))
  (uuid:uuid-to-integer u))

(defmethod integer-of ((s string))
  (integer-of (vector-of s)))

;; ------------

(defmethod vector-of ((v vector))
  v)

(deftype vector-ub8 ()
  `(simple-array (unsigned-byte 8) (*)))

(defun string-to-vector (str)
  (map 'vector-ub8 'char-code str))

(defmethod vector-of ((s string))
  (string-to-vector s))

(defmethod vector-of ((u uuid:uuid))
  (uuid:uuid-to-byte-array u))

(defmethod vector-of ((x integer))
  (convert-integer-to-vector x))

;; -----------------------------------------------

(defun make-nbit-prime (nbits hash-type)
  (assert (>= (* 8 (ironclad:digest-length hash-type)) nbits))
  (let* ((lower (ash 1 nbits))
         (upper (+ lower lower))
         (mask  (1- (ash lower -1)))
         (base  (+ mask 2)))
    (loop for seed = (random-between lower upper)
          for u = (let ((dig (ironclad:make-digest hash-type)))
                    (ironclad:update-digest dig (vector-of seed))
                    (logand (integer-of (ironclad:produce-digest dig)) mask))
          for q = (- (+ base u) (logand u 1))
          when (is-prime? q)
          do (return (values q seed))) ))

(defun gen-dsa-primes (nbits lbits hash-type)
  (let* ((hbits (* 8 (ironclad:digest-length hash-type)))
         (hfact (ash 1 hbits))
         (n     (1- (ceiling lbits hbits)))
         (b     (- lbits 1 (* n hbits)))
         (bmask (1- (ash 1 b)))
         (2^lm1 (ash 1 (1- lbits)))
         (hmask (1- (ash 1 nbits))))
    (loop for (q seed) = (multiple-value-list (make-nbit-prime nbits hash-type))
          do
          (format t "~&Q = ~A" q)
          (loop for counter from 0 below (* 4 lbits)
                for offset from 1 by (1+ n)
                for v = (make-array (1+ n))
                do
                (loop for j from 0 to n do
                      (setf (aref v j) (let ((dig (ironclad:make-digest hash-type)))
                                         (ironclad:update-digest dig
                                                                 (vector-of
                                                                  (logand (+ seed offset j)
                                                                          hmask)))
                                         (integer-of (ironclad:produce-digest dig)) )))
                (let* ((w (loop for ix from 0 to n
                                for k = 1 then (* k hfact)
                                sum (* k (let ((x (aref v ix)))
                                           (if (= ix n)
                                               (logand x bmask)
                                             x)))))
                       (x (+ w 2^lm1))
                       (c (mod x (* 2 q)))
                       (p (- x (1- c))))
                  (when (and (>= p 2^lm1)
                             (is-prime? p))
                    (return-from gen-dsa-primes (values p q)) ))) )))
                  

;; -------------------------------------------------------

(defun make-2kp+1-prime (nbits &optional (mr-iters 50) (kdiff 10))
  ;; a safe prime is of the form (2*k*U+1) where U is a large prime, and k is a small multiplier.
  ;; NLEVELS indicates how many component (P-1) should be strong
  ;;  NLEVELS: T   = all levels should be strong and proven prime
  ;;           NIL = none should be strong -> probabilistic prime
  ;;           <integer> = nbr of strong component levels - proven or probabilistic result
  ;;                 If the printout shows: Make Base Prime: Absolute
  ;;                    then we will have a proven prime result.
  ;;                 Otherwise, we will have a probabilistic prime result.
  (let ((base   (ash 1 (- nbits kdiff)))
        (kstart (ash 1 (- (max kdiff 2) 2))))
    (loop for u = (make-prime base mr-iters) do
          (let ((incr (+ u u)))
            (um:nlet-tail iter ((p (1+ (* kstart incr))))
              (let ((nb (integer-length p)))
                (cond ((< nb nbits) (iter (+ p incr)))
                      ((= nb nbits) (if (is-prime? p mr-iters)
                                        (return-from make-2kp+1-prime p)
                                      (iter (+ p incr))))
                      (t nil)) )) )) ))

;; ------------------------------------------------------------------------
;;
#|
;; 128-bit primes of the form 2*q+1, q prime

110733130861390325273223714915708815940841239025086956492482512302431841768612420553926491015864569438667372511804512595013714067655948701888623803535130698355087737510353066609832316632241419637710746726030623432724836837829706946251221511561676095077906242421212235304835973028416739591897787588131953372219

110733130861390325273223714915708815940841239025086956492482512302431841768612420553926491015864569438667372511804512595013714067655948701888623803535130698355087737510353066609832316632241419637710746726030623432724836837829706946251221511561676095077906242421212235304835973028416739591897787588131953542247

110733130861390325273223714915708815940841239025086956492482512302431841768612420553926491015864569438667372511804512595013714067655948701888623803535130698355087737510353066609832316632241419637710746726030623432724836837829706946251221511561676095077906242421212235304835973028416739591897787588131954440903

110733130861390325273223714915708815940841239025086956492482512302431841768612420553926491015864569438667372511804512595013714067655948701888623803535130698355087737510353066609832316632241419637710746726030623432724836837829706946251221511561676095077906242421212235304835973028416739591897787588131954912347

110733130861390325273223714915708815940841239025086956492482512302431841768612420553926491015864569438667372511804512595013714067655948701888623803535130698355087737510353066609832316632241419637710746726030623432724836837829706946251221511561676095077906242421212235304835973028416739591897787588131955850507

121636686180494859351257062932794606637412604912349654946895786579232564207004269642515857596415775103556609593312343761619399270973617792375942737442895269169033370326978057626223123112060577736588615703266513961275805218656107856015364524626128664217938276131801578811416035728088503132610052344318901362059

121636686180494859351257062932794606637412604912349654946895786579232564207004269642515857596415775103556609593312343761619399270973617792375942737442895269169033370326978057626223123112060577736588615703266513961275805218656107856015364524626128664217938276131801578811416035728088503132610052344318901390307

121636686180494859351257062932794606637412604912349654946895786579232564207004269642515857596415775103556609593312343761619399270973617792375942737442895269169033370326978057626223123112060577736588615703266513961275805218656107856015364524626128664217938276131801578811416035728088503132610052344318901683443

121636686180494859351257062932794606637412604912349654946895786579232564207004269642515857596415775103556609593312343761619399270973617792375942737442895269169033370326978057626223123112060577736588615703266513961275805218656107856015364524626128664217938276131801578811416035728088503132610052344318903543947

121636686180494859351257062932794606637412604912349654946895786579232564207004269642515857596415775103556609593312343761619399270973617792375942737442895269169033370326978057626223123112060577736588615703266513961275805218656107856015364524626128664217938276131801578811416035728088503132610052344318903574907

121636686180494859351257062932794606637412604912349654946895786579232564207004269642515857596415775103556609593312343761619399270973617792375942737442895269169033370326978057626223123112060577736588615703266513961275805218656107856015364524626128664217938276131801578811416035728088503132610052344318904076759

125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701859684739

125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701859958483

125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701860004983

125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701860519219

125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701860543159

125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701860928467
 
125874072370696794236809845433085663986098930803397897971969785329075699829719788543438626417022264395405304793095682096307356853890121801628891244488742033671712430247054683916845968149901694358504515381754379369711976936208218925608516005824179148210832986948494616110927217460975273585946238026701862024703

126109499319730975072643891561817777090078234159970968064380371675130087731049251176894058682887645606400141915814997634886187326048858326430138948012633141618386065743242120967800116033859745494421234716111317867182827099082997052278782304518545867425519912626164985878601696040793159811247884048409248818959

128491880968628397295439786826090448071799216337734462740114772282882719015462829336458309749315736156755660134536485137110659061247296862257126145128366477583072943793395168181329513531307158492226048781300219481587675392558457213959505574314676368201921410594701570918097602501976047163949823423891177344307

141748196247699632500986808906971058611550869395507694922590511654883610826606305709233518138170423599408401107955728398351083656546276901782102104779437208477121080540723797281574158246955100480326947890864922029402193565967828467664068413921533244861976824769684937384113548662768837624008184082755291790299

141748196247699632500986808906971058611550869395507694922590511654883610826606305709233518138170423599408401107955728398351083656546276901782102104779437208477121080540723797281574158246955100480326947890864922029402193565967828467664068413921533244861976824769684937384113548662768837624008184082755292134987

141748196247699632500986808906971058611550869395507694922590511654883610826606305709233518138170423599408401107955728398351083656546276901782102104779437208477121080540723797281574158246955100480326947890864922029402193565967828467664068413921533244861976824769684937384113548662768837624008184082755293540883

141748196247699632500986808906971058611550869395507694922590511654883610826606305709233518138170423599408401107955728398351083656546276901782102104779437208477121080540723797281574158246955100480326947890864922029402193565967828467664068413921533244861976824769684937384113548662768837624008184082755293652723

141748196247699632500986808906971058611550869395507694922590511654883610826606305709233518138170423599408401107955728398351083656546276901782102104779437208477121080540723797281574158246955100480326947890864922029402193565967828467664068413921533244861976824769684937384113548662768837624008184082755294155379

148757362201947503107350056672229226015846526670066439700592774879557900230840629657269570465529206476597033680891372798894319849683448592503697981688782652297297277079736508164805831459847318748821103794744570332163018109903223211428324274319113375615993139362231079938382046507249340896843828581789856351283
 
148757362201947503107350056672229226015846526670066439700592774879557900230840629657269570465529206476597033680891372798894319849683448592503697981688782652297297277079736508164805831459847318748821103794744570332163018109903223211428324274319113375615993139362231079938382046507249340896843828581789856714019

148757362201947503107350056672229226015846526670066439700592774879557900230840629657269570465529206476597033680891372798894319849683448592503697981688782652297297277079736508164805831459847318748821103794744570332163018109903223211428324274319113375615993139362231079938382046507249340896843828581789859728239

148757362201947503107350056672229226015846526670066439700592774879557900230840629657269570465529206476597033680891372798894319849683448592503697981688782652297297277079736508164805831459847318748821103794744570332163018109903223211428324274319113375615993139362231079938382046507249340896843828581789859801019

171050948217267985673436464661659295399250022438505679154384395855169504304692006719499500515962744962166121514204840502682644185949107664427757974606706183917926195232484819082843189438596224733882477303185478371183408054236233026031072171840648622216400655791611595902476601403117316676079943273588450078823

171050948217267985673436464661659295399250022438505679154384395855169504304692006719499500515962744962166121514204840502682644185949107664427757974606706183917926195232484819082843189438596224733882477303185478371183408054236233026031072171840648622216400655791611595902476601403117316676079943273588450196183

171050948217267985673436464661659295399250022438505679154384395855169504304692006719499500515962744962166121514204840502682644185949107664427757974606706183917926195232484819082843189438596224733882477303185478371183408054236233026031072171840648622216400655791611595902476601403117316676079943273588451879867

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645108962583

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645111589323

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645111835839

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645112131423

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645112476843

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645112921323

172995140565004716223251173179587368013743810348742732185154436072499697094058486203981224807130120428602455987779848226007562591515631828603480639245218673811592117300197954957487164919194991946812362886425360148736764108748693333982655118403291835063746524836780814704300732233146047676283200064645112941687

94151763009122605199016046651212124826776428399875868337139174973411137492045166520124166776862436148366239370758238172290409368422703069135426261530234619164046611105339201384911749696986955053718391402729253242954270262493136526039898838280505289291728617563351299061826021401294993431392010748869152068839

;; ------------------------------------------
;; 3072-bit primes
3831312298414591100618585476831656081302761106445333367964445766101605673062226600081798851596800064959373358447986085212005152872491415633440971408310613696416816247131251191127605769187849632722625692595830232912120509411459227561219422278588388521279190819787228640437741365583664399535445890926433594017171015279994876446437819329266951171495721820598600803760517032910294611120426233980568547605244805692942542173452077092113363027552652352478710408874217449160078586013212918771714470031143005776430786677042471224235686920661387334363033446494987716085736076096310420090635775820390682772337164190737858010311492171925311511494388679285585265392735735001153369804939762264859453754474981785662042445364307068663290064088038175751118437060927313568593020052845029608309540759305204515401915036636381115496655888777888350958428841118551928804350704193370950059587452549447387500168231728657073005226795752366621721650043
 
3831312298414591100618585476831656081302761106445333367964445766101605673062226600081798851596800064959373358447986085212005152872491415633440971408310613696416816247131251191127605769187849632722625692595830232912120509411459227561219422278588388521279190819787228640437741365583664399535445890926433594017171015279994876446437819329266951171495721820598600803760517032910294611120426233980568547605244805692942542173452077092113363027552652352478710408874217449160078586013212918771714470031143005776430786677042471224235686920661387334363033446494987716085736076096310420090635775820390682772337164190737858010311492171925311511494388679285585265392735735001153369804939762264859453754474981785662042445364307068663290064088038175751118437060927313568593020052845029608309540759305204515401915036636381115496655888777888350958428841118551928804350704193370950059587452549447387500168231728657073005226795752366621722201983


4958785161367786818650506991217893560732160987155560166147355256318343730672007000617972370099971962494700856552987951937218247097059454941629097569639188036827211147105738076734511317020078487276286491358372632084332882695850067691398348189906787233169219752059424030467571339174639350733762588153818101491313128431544377463490330085235988821840800868743767765574584051564101839230517392257253335131309484284694143823340149541779335482773415887622206258688531561916037338908636837083115429203108131090935827916738395997365880322515121312964782143228810034285376087076317654664015150046265987030969619227677703418063968963098558420966734477221970070843218349274881138321331257168869525821381696480724905193809829613249838623332872950546090941864063816507225702124570173612604409545843018398456301621339452868044785522603462621441423561331069806338876113287105372571110488493829754039066151606455764564313406778878957833765827

|#

(defun make-2q+1-prime (kbits &optional (mr-iters 50))
  ;; If P = 2*Q+1, for Q prime, and all higher primes are 6*k+/-1, then
  ;;   if Q = 6*k-1 then P = 12*k-1
  ;;   if Q = 6*k+1 then P = 12*k+3, but that is divisible by 3, and so cannot be.
  ;; Therefore we want to look only for Q = 6*k-1, Q prime, and P = 2*Q+1 prime
  #F
  (declare (integer kbits))
  (let* ((lower (1+ (ash 1 (- kbits 2))))
         (upper (1- (ash 1 (1- kbits))))
         (niter 0)
         (range kbits))
    (declare (integer lower upper)
             (fixnum niter range))
    
    ;; Since fundamental theorem states that pi(n) ~ n/ln(n),
    ;; or density of primes is inversely proportional to number of decimal digits,
    ;; on average about ln(k)/2 trials will be needed

    ;; If P = 2*Q+1, for Q prime, and all higher primes are 6*k+/-1, then
    ;;   if P = 6*k-1 = 2*Q+1 ==> Q = 3*k-1, and Q prime
    ;;   if P = 6*k+1 = 2*Q+1 ==> Q = 3*k, so Q can't be prime
    ;; Therefore, examine only the candidates q = 3*k-1.

    (labels ((test (n)
               (declare (integer n))
               (and (is-prime? n mr-iters)
                    (let ((p (1+ (ash n 1))))
                      (and (is-prime? p mr-iters)
                           p)))))

      (loop for p =
            (let* ((x  (ash 1 (- kbits 2)))
                   (y  (logior x (mt-random x)))
                   (r  (mod (- 5 y) 6))
                   (n  (+ y r)))
              (declare (integer n x y)
                       (fixnum r))
              (print n)
              (or
               (um:nlet-tail iter> ((n  n)
                                    (ix 0))
                 (declare (integer n)
                          (fixnum ix))
                 (when (< ix range)
                   (incf niter)
                   (when (< n upper)
                     (um:if-let (p (test n))
                         p
                       (iter> (+ n 6) (1+ ix))) )))

               (um:nlet-tail iter< ((n (- n 6))
                                    (ix 0))
                 (declare (integer n)
                          (fixnum ix))
                 (when (< ix range)
                   (incf niter)
                   (when (> n lower)
                     (um:if-let (p (test n))
                         p
                       (iter< (- n 6) (1+ ix)) ))))))
              when p do (format t "~&niter = ~A" niter)
              when p return p) )))

;; -------------------------------------------------------------------

(defun sqr (n)
  (* n n))

(defun big-log (n)
  ;; log for BigNums
  (* (integer-length n) (log 2)))

(defun pdens (n)
  ;; approx local density of primes near N
  ;; from Bertrand's postulate
  (/ (big-log n)))


(defun sieve-2q+1 (nbits &optional (mr-iters 50))
  ;; If P = 2*Q+1, for Q prime, and all higher primes are 6*k+/-1, then
  ;;   if Q = 6*k-1 then P = 12*k-1
  ;;   if Q = 6*k+1 then P = 12*k+3, but that is divisible by 3, and so cannot be.
  ;; Therefore we want to look only for Q = 6*k-1, Q prime, and P = 2*Q+1 prime
  #F
  (declare (fixnum nbits mr-iters))
  (let* ((tbits 18) ;; big enough to search 1.5M candidates
         (ns    (ash 1 tbits))
         (arr   (make-array ns
                            :element-type    'bit
                            :initial-element 0))
         (base  (logior (ash 1 (- nbits 2)) ;; MSB
                        (ash (mt-random (ash 1 (- nbits 3 tbits 3))) (+ tbits 3)) ;; body
                        (ash (mt-random 3) tbits))) ;; 2 remainder of 8 by 6
         (r     (mod (- 5 base) 6))
         (lower (+ base r))) ;; 2^(nbits-1) < 2*lower+1 < 2*lower+1+12*ns < 2^nbits
         
    ;; We have lower = -1 mod 6 = 5 mod 6, or lower = 6*m-1
    ;; each bit in arr represents (lower + 6*k) = (6*(k+m)-1)
    ;; First mark the ones that are known composite numbers from our *primes* table.
    (declare (fixnum ns kmod3)
             (integer vmin lower)
             ((array bit (*)) arr))

    (dolist (p (cddr *primes*))
      (declare (fixnum p))
      (do ((offs (the fixnum
                      (let ((h (mod 6 p)))
                        (declare (fixnum h))
                        (mod (- (* (mod lower p) (inv-mod h p))) p)))
                 (the fixnum (+ offs p))))
          ((>= offs ns))
        (declare (fixnum offs))
        (setf (aref arr offs) 1)))

    ;; Now scan for possible primes -- not marked composite
    ;; Form q = (lower + 6*k) and test for q prime, then test for p = 2*q+1 prime
    (do ((k  0     (1+ k))
         (ps nil))
        ((>= k ns) (nreverse ps))
      (declare (fixnum k))
      (when (zerop (aref arr k))
        (let ((q (+ lower (* 6 k))))
          (declare (integer q))
          (when (is-prime? q mr-iters)
            (let ((p (1+ (ash q 1))))
              (declare (integer p))
              (when (is-prime? p mr-iters)
                (push p ps)) ))))) ))

