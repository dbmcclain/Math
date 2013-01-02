
;; ------------------------------------------------------------
;; Root finding via Newton's Method
;;
(defun abserror (eps)
  (lambda (a b)
    (<= (abs (- a b)) eps)))

(defun relerror (eps)
  (let ((aerr (abserror eps)))
    (lambda (a b)
      (if (zerop b)
          (<= (abs a) eps)
        (funcall aerr (/ a b) 1.0d0)))
    ))

(defun zbrac (fn xmin xmax
                 &key (maxdepth 4))
  (labels ((iter (niter x1 f1 x2 f2)
             (when (< niter maxdepth)
               (cond ((zerop f1) (list x1 x1))
                     ((zerop f2) (list x2 x2))
                     ((not (plusp (* f1 f2))) (list x1 x2))
                     (t (let* ((xmid  (* 0.5d0 (+ x1 x2)))
                               (fxmid (funcall fn xmid)))
                          (or (iter (1+ niter) x1 f1 xmid fxmid)
                              (iter (1+ niter) xmid fxmid x2 f2))
                          ))
                     ))
             ))
    (iter 0
          xmin (funcall fn xmin)
          xmax (funcall fn xmax))
    ))

(defparameter *default-convergence-test* (abserror 1.0d-8))

(defun newton-root (fn derfn xmin xmax
                       &key
                       (test     *default-convergence-test*)
                       (maxiters 20) ;max nbr of iterations permitted
                       (maxdepth  4) ;max depth for zbrac root bracketing
                       (eps      1.0d-8) ;function value at root must be within this of zero
                       start ;starting guess for root
                       &aux
                       (x0  (or start
                                (um:if-let (pair (zbrac fn xmin xmax :maxdepth maxdepth))
                                           ;; choose the endpoint with the steepest slope
                                           (destructuring-bind (x1 x2) pair
                                             (let ((der1 (abs (funcall derfn x1)))
                                                   (der2 (abs (funcall derfn x2))))
                                               (if (>= der1 der2)
                                                   x1
                                                 x2))
                                             ))))
                       (xmn (min xmin xmax))
                       (xmx (max xmin xmax)))

  ;; Return (values NIL)         if we couldn't bracket a root in the interval [xmin, xmax]
  ;; Return (values NIL niter)   if we can't converge within maxiters iterations
  ;; Return (values xroot niter) if we succeed
  
  (labels ((sharpen (x xprev)
             (let ((f   (funcall fn    x))
                   (der (funcall derfn x)))
               
               (cond ((zerop der)
                      ;; can't divide by zero derivative
                      ;; try again using midpoint of x and xprev
                      (* 0.5d0 (+ x xprev)))
                     
                     (t
                      ;; the usual case
                      ;; perform one round of Newton sharpening
                      ;; but restrict result to interval [xmin, xmax]
                      (min xmx (max xmn (- x (/ f der)))))
                     ))
             ))

    (when x0
      ;; else couldn't bracket root in zbrac
      (loop for niter from 0
            for (x xprev) = (list x0 nil) then (list (sharpen x xprev) x)
            do
            (cond (;; successful exit
                   (and xprev
                        (funcall test x xprev)
                        (< (abs (funcall fn x)) eps))
                   (return (values x niter)))
                  
                  (;; nonconvergence failure
                   (> niter maxiters)
                   (return (values nil niter)))
                  )))
    ))

;; ----------------------------------------------------------------------------
;;
(defun zbrac-outward (fn x1 x2)
  (let* ((f1 (funcall fn x1))
         (f2 (funcall fn x2))
         (logx1 (log x1))
         (logx2 (log x2)))
    (loop for ix from 1
          until (minusp (* f1 f2)) do
          ;;(print (list (exp logx1) f1 (exp logx2) f2))
          (if (> ix 50)
              (error "Bad initial range in zbrac-outward"))
          (if (< (abs f1) (abs f2))
              (setf f1 (funcall fn (exp (incf logx1 (* 1.6 (- logx1 logx2))))))
            (setf f2 (funcall fn (exp (incf logx2 (* 1.6 (- logx2 logx1)))))))
          )
    (values (exp logx1) f1 (exp logx2) f2)))
          
(defun ridder (fn xlow xhi &optional (epsilon 1.0e-6))
  ;; Ridder's superlinear method for finding roots
  (multiple-value-bind (xlow flow xhi fhi)
      (zbrac-outward fn xlow xhi)
    (if (zerop flow)
        xlow
      (if (zerop fhi)
          xhi
        (let* ((lxlow (log xlow))
               (lxhi  (log xhi))
               (ans  1.11e30))
          (labels
              ((ridder-iter (niter)
                 (when (> niter 60)
                   (error "Ridder: no root found after 60 iterations"))
                 (let* ((lxmid (* 0.5d0 (+ lxlow lxhi)))
                        (fmid (funcall fn (exp lxmid)))
                        (s  (sqrt (- (* fmid fmid) (* flow fhi)))))
                   (if (zerop s)
                       (exp ans)
                     (let* ((lxnew (+ lxmid
                                      (* (- lxmid lxlow)
                                         (let ((ratio (/ fmid s)))
                                           (if (>= flow fhi)
                                               ratio
                                             (- ratio)))
                                         ))))
                       (if (<= (abs (- lxnew (shiftf ans lxnew))) epsilon)
                           (exp lxnew)
                         (let ((fnew (funcall fn (exp lxnew))))
                           (if (zerop fnew)
                               (exp lxnew)
                             (progn
                               ;;(print (list xnew fnew))
                               (cond ((minusp (* fnew fmid))
                                      (setf lxlow lxmid
                                            lxhi  lxnew
                                            flow fmid
                                            fhi  fnew))
                                     ((minusp (* flow fnew))
                                      (setf lxhi lxnew
                                            fhi fnew))
                                     ((minusp (* fhi fnew))
                                      (setf lxlow lxnew
                                            flow fnew))
                                     (t (error "Ridder: Can't happen")))
                               (if (<= (abs (- lxhi lxlow)) epsilon)
                                   (exp lxnew)
                                 (ridder-iter (1+ niter)))))
                           ))))
                   )))
            (ridder-iter 1)))
        ))))

(defun root-by-bisection (fn x0 x1 &optional (epsilon 1.0d-10))
  (multiple-value-bind (x1 y1 x2 y2)
      (zbrac-outward fn x0 x1)
    (declare (ignore y2))
    (let ((logx1 (log (max 0.01 x1)))
          (logx2 (log (max 0.01 x2))))
      (loop while (> (abs (- logx1 logx2)) epsilon) do
            (let* ((xmid (exp (* 0.5d0 (+ logx1 logx2))))
                   (ymid (funcall fn xmid)))
              ;;(print (list xmid fmid))
              (if (not (plusp (* y1 ymid)))
                  (setf logx2 (log xmid))
                (setf logx1 (log xmid)
                      y1 ymid))
              ))
      (exp logx1))))

(defun find-root (fn xlow xhi &optional (epsilon 1.0d-10))
  ;;(ridder fn xlow xhi epsilon)
  (root-by-bisection fn xlow xhi epsilon)
  ;;(newton-root fn (nderfn fn) xlow xhi)
  )

;; ----------------------------------------------------------------
#|
(defun fnrise (n)
  (labels ((fn (n)
             (* 0.06 (expt 0.94 n))))
    (loop for ix from 0 to n sum (fn ix))))

(let ((win (plt:wset 0)))
  (plt:clear win)
  (plt:plot win (loop for ix from 0 to 100 collect (fnrise ix)))
  (plt:plot win (loop for ix from 0 to 100 collect (- 1.0 (expt 0.94 ix)))))

(find-root (lambda (alpha) (- 0.99 (- 1.0 (expt (- 1.0 alpha) 193)))) 1e-6 1)
     

|#
