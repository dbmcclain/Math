;; Remez.lisp -- Remez's method for finding minimax polynomial and rational approximations
;; DM/RAL 01/07
;;

#|
(defun quad-find-peak (ym1 y0 yp1)
  ;; for equispaced data
  ;; quadratic interpolation to locate the min or max
  (/ (* 0.5 (- ym1 yp1)) (- (+ ym1 yp1) (+ y0 y0))))
|#

(defun abserror (a b)
  (abs (- a b)))

(defun relerror (a b)
  (if (zerop b)
      (abs a)
    (abserror (/ a b) 1.0d0)))

#|
;; these can produce instabilities
(defun quad-find-peak (x0 x1 x2 y0 y1 y2)
  (let* ((d0 (/ y0 (* (- x0 x1) (- x0 x2))))
         (d1 (/ y1 (* (- x1 x0) (- x1 x2))))
         (d2 (/ y2 (* (- x2 x0) (- x2 x1))))
         (n0 (* d0 (+ x1 x2)))
         (n1 (* d1 (+ x0 x2)))
         (n2 (* d2 (+ x0 x1))))
    (/ (+ n0 n1 n2) (+ d0 d1 d2) 2.0d0)
    ))

(defun find-extremum (fn a b c fa fb fc)
  ;; using iterated quadratic interpolation
  ;; on entry (a < b < c) and (fa < fb) and (fc < fb)
  (um:with-tail-pure-code
    (assert (< a b c))
    (assert (<= fa fb))
    (assert (<= fc fb))
    (let* ((m   (quad-find-peak a b c fa fb fc))
           (fm  (funcall fn m)))
      ;;(print (list a m1 b m2 c))
      (cond
       ((< (relerror fm fb) 1.0d-8) m)
       ((and (< a m b)
             (< fa fm)
             (< fb fm))  (find-extremum fn a m b fa fm fb))
       ((and (< b m c)
             (< fb fm)
             (< fc fm))  (find-extremum fn b m c fb fm fc))
       (t (error "Can't happen in FIND-EXTREMUM ~A" (list a b c m)))
       ))))
|#

#| |#
(defun find-extremum (fn a b c fb)
  ;; using interval bisection
  ;; a < b < c && fn(a) < fn(b) && fn(c) < fn(b)
  (um:with-tail-pure-code
    (labels ((iter (a b c fb)
               (let* ((m1  (* 0.5d0 (+ a b)))
                      (fm1 (abs (funcall fn m1))))
                 (cond
                  ((< (abserror m1 b) 1.0d-8) m1)
                  ((< (relerror fm1 fb) 1.0d-8) m1)
                  ((> fm1 fb) (iter a m1 b fm1))
                  (t  (let* ((m2  (* 0.5d0 (+ b c)))
                             (fm2 (abs (funcall fn m2))))
                        (cond ((< (abserror m2 b) 1.0d-8) m2)
                              ((< (relerror fm2 fb) 1.0d-8) m2)
                              ((> fm2 fb) (iter b m2 c fm2))
                              ((and (> fb fm1)
                                    (> fb fm2)) (iter m1 b m2 fb))
                              (t (error "Can't happen in FIND-EXTREMUM"))
                              )))
                  ))))
      (iter a b c (abs fb))
      )))
#| |#

(defun find-nodes (fn xs ys)
  (let* ((limit (1- (length ys)))
         (xsv   (coerce xs 'vector))
         (ysv   (coerce ys 'vector))
         (prev  (aref ysv 0))
         (nodes nil))
    (loop for ix from 1 below limit do
          (let* ((ym1 (aref ysv (1- ix)))
                 (y0  (aref ysv ix))
                 (yp1 (aref ysv (1+ ix))))
            ;; potential extremum on opposite side of zero from previous?
            (when (and (minusp (* y0 prev))
                       (>= (abs y0) (abs ym1))
                       (>= (abs y0) (abs yp1)))
              (push (find-extremum fn
                                   (aref xsv (1- ix))
                                   (aref xsv ix)
                                   (aref xsv (1+ ix))
                                   y0)
                    nodes)
              (setf prev y0))
            ))
    ;; list of new nodes sandwitched between first and last abscissae
    (cons (first xs)
          (append (nreverse nodes)
                  (last xs)))))

(defun solve-for-coffs (nodes ys wts np nq)
  ;; solve the pseudo-linear system for the coefficients
  (um:with-tail-pure-code
    (let ((ma    (make-array (list (+ np nq 2)
                                   (+ np nq 2))
                             :element-type 'double-float)))
      (labels ((iter (eps coffs niter)
                 (if (> niter 20)
                     (subseq coffs 0 (+ np nq 1))
                   (progn
                     (loop for jx from 0 to (+ np nq 1)
                           and xjx in nodes
                           and yjx in ys
                           and wjx in wts
                           do
                           (loop for kx from 0 to np do
                                 (setf (aref ma jx kx) (expt xjx kx)))
                           (loop for kx from 1 to nq do
                                 (setf (aref ma jx (+ np kx))
                                       (- (* (expt xjx kx)
                                             (+ yjx
                                                (/ (* (expt -1 jx) eps) wjx))
                                             ))))
                           (setf (aref ma jx (+ np nq 1))
                                 (- (/ (expt -1 jx) wjx))))
                     (let* ((coffs (vm:dgesvd-solve ma
                                                    (coerce ys 'vector)
                                                    :tolerance 0))
                            (eps-new (aref coffs (+ np nq 1))))
                       (if (< (relerror eps eps-new) 1.0d-4)
                           (subseq coffs 0 (+ np nq 1))
                         (iter eps-new coffs (1+ niter))
                         ))))))
        (iter 0.0d0 nil 1))
      )))

(define-condition bad-error-curve (error) ()
  (:report (lambda (condition stream)
             (declare (ignore condition))
             (princ "Nonstandard Error Curve in REMEZ" stream))
   ))

(defun weighted-minimax-rational-approx (fn domain fnwt order)
  (um:with-tail-pure-code
    (destructuring-bind (np nq) order
      (destructuring-bind (xmin xmax) domain
        (let* ((rrfn (lambda (x)
                       (/ (- (* 2d0 x) (+ xmax xmin))
                          (- xmax xmin))))
               (irrfn (lambda (x)
                        (* 0.5d0 (+ (* x (- xmax xmin))
                                    (+ xmax xmin)))))
               (xs    (loop for ix from -100 to 100 collect
                            (* 0.01d0 ix)))
               (xxs   (mapcar irrfn xs))
               (ys    (mapcar fn xxs))
               (wts   (mapcar fnwt xxs)))

          (labels ((solve-for-initial-coffs (nodes)
                     (let* ((xns (mapcar irrfn nodes))
                            (yns (mapcar fn xns))
                            (wtns (mapcar fnwt xns)))
                       (solve-for-coffs nodes yns wtns np nq)))
                   
                   (eval-coffs (coffs x)
                     (let* ((num (loop for ix from 0 to np sum
                                       (* (aref coffs ix)
                                          (expt x ix))))
                            (den (loop for ix from 1 to nq sum
                                       (* (aref coffs (+ ix np))
                                          (expt x ix)))))
                       (/ num (+ 1.0d0 den))))
                                  
                   (fnerr (coffs)
                     (lambda (x)
                       (let ((xx (funcall irrfn x)))
                         (* (funcall fnwt xx) (- (eval-coffs coffs x) (funcall fn xx)))
                         )))
                   
                   (check-coffs (nodes coffs)
                     (let* ((yhats (mapcar (lambda (x)
                                             (eval-coffs coffs x))
                                           xs))
                            (errs  (mapcar #'- yhats ys))
                            (errsw (mapcar #'* errs wts))
                            
                            (xns   (mapcar irrfn nodes))
                            (yns   (mapcar fn xns))
                            (wtns  (mapcar fnwt xns))
                            (errns (mapcar (lambda (x y)
                                             (- (eval-coffs coffs x) y))
                                           nodes yns))
                            (errnsw (mapcar #'* errns wtns)))
                       
                       (let ((win (plt:wset 1 :clear t)))
                         (plt:plot win xs errsw)
                         (plt:plot win nodes errnsw
                                   :symbol :circle
                                   :color  :blue)
                         
                         (let* ((nodes-new (find-nodes (fnerr coffs) xs errsw))
                                (xns-new   (mapcar irrfn nodes-new))
                                (yns-new   (mapcar fn xns-new))
                                (wtns-new  (mapcar fnwt xns-new))
                                (errns-new (mapcar (lambda (x y)
                                                     (- (eval-coffs coffs x) y))
                                                   nodes-new yns-new))
                                (errnsw-new (mapcar #'* errns-new wtns-new)))
                           
                           (plt:plot win nodes-new errnsw-new
                                     :symbol :square
                                     :color  :red)
                           
                         (let* ((is-stderr (= (length nodes-new) (+ np nq 2)))
                                (eabs      (mapcar #'abs errnsw-new))
                                (emax      (reduce #'max eabs))
                                (emin      (reduce #'min eabs)))
                           
                           (values nodes-new is-stderr emax (- emax emin)))
                         ))))

                   (iter-remez (nodes coffs niter)
                     (multiple-value-bind (nodes-new is-stderr err derr)
                         (check-coffs nodes coffs)
                       (unless is-stderr
                         (error 'bad-error-curve))
                       ;;(sleep 0.5)
                       (let ((coffs-new (let* ((xns  (mapcar irrfn nodes-new))
                                               (yns  (mapcar fn xns))
                                               (wtns (mapcar fnwt xns)))
                                          (solve-for-coffs nodes-new yns wtns np nq))))
                         (if (> niter 20)
                             ;;(every (lambda (a b)
                             ;;         (< (abserror a b) (* 0.1d0 err)))
                             ;;       coffs coffs-new)
                             (list :nodes nodes-new
                                   :coffs coffs-new
                                   :order order
                                   :err   err
                                   :derr  (/ derr err)
                                   :rrfn  rrfn)
                           (iter-remez nodes-new coffs-new (1+ niter)))
                         )))
                   
                   (iter (nodes niter)
                     (let ((coffs (solve-for-initial-coffs nodes)))
                       (handler-case
                           (iter-remez nodes coffs 1)
                         (bad-error-curve ()
                           (if (> niter 50)
                               (error 'bad-error-curve)
                             (let ((new-nodes (loop for node in nodes
                                                    and ix from 0
                                                    collect
                                                    (if (or (zerop ix)
                                                            (= ix (+ np nq 1)))
                                                        node
                                                      (+ node (random 0.1d0) -0.05d0))
                                                    )))
                               (print "retry with dithered nodes")
                               (iter new-nodes (1+ niter)))
                             )))
                       )))

            (let ((nodes (loop for ix from 0 to (+ np nq 1) collect
                               (cos (/ (* pi ix) (+ np nq 1))))))
              (iter nodes 1))
            ))
        ))
    ))

(defun minimax-rational-approx (fn domain order)
  ;; based on absolute error
  (weighted-minimax-rational-approx fn domain
                                    (lambda (x)
                                      (declare (ignore x))
                                      1.0d0)
                                    order))

(defun relative-minimax-rational-approx (fn domain order)
  ;; based on relative error
  (weighted-minimax-rational-approx fn domain
                                    (lambda (x)
                                      (/ (funcall fn x)))
                                    order))


#|
(minimax-rational-approx #'exp '(0.0d0 1.0d0) '(4 0))
(minimax-rational-approx (lambda (x) (log x 2)) '(0.5d0 1.0d0) '(2 0))
(minimax-rational-approx (lambda (x) (expt 2 x)) '(0.0d0 1.0d0) '(2 0))
(minimax-rational-approx (lambda (x) (expt 2 x)) '(-1.0d0 1.0d0) '(3 0))
(minimax-rational-approx #'exp '(0.0d0 1.0d0) '(4 4))

|#


