
(proclaim '(debug 2)) ;; eliminate tail recursion

(defun my-iexpt (x n)
  (labels ((iter (x n ans)
             (if (zerop n)
                 ans
               (progn
                 ;; (format t "~%x=~F, n=~D, ans=~F"
                 ;;         (float x 1.0d0) n (float ans 1.0d0))
                 (iter (* x x) (ash n -1)
                       (if (zerop (logand n 1))
                           ans
                         (* ans x))))
               )))
    (iter x n 1)
    ))

(defun t4x ()
  (let* ((x (loop for ix from 1 to 128 collect
                  (+ 1.0F0 (/ (float ix 1.0F0) 65536.0F0))))
         (n (loop for ix from 1 to 128 collect
                  (+ 65500 (* 2 ix))))
         (x3n (loop for x in x and n in n collect
                    (expt x (* 3 n))))
         (xxxn (loop for x in x and n in n collect
                     (expt (* x x x) n)))
         (actualx (loop for ix from 1 to 128 collect
                        (+ 1.0d0 (/ (float ix 1.0d0) 65536.0d0))))
         (actualx3n (loop for x in actualx and n in n collect
                          (expt x (* 3 n))))
         (adiffx3n (/ (vmax (mapcar #'(lambda (x3n ref)
                                        (abs (/ (- x3n ref) ref)))
                                    x3n actualx3n))
                      eps))
         (adiffxxxn (/ (vmax (mapcar #'(lambda (xxxn ref)
                                         (abs (/ (- xxxn ref) ref)))
                                     xxxn actualx3n))
                       eps))
         (diff (/ (vmax (mapcar #'(lambda (x3n xxxn)
                                    (abs (/ (- x3n xxxn) x3n)))
                                x3n xxxn))
                  eps)))
    ;;(break)
    (pnl (format nil "y^N errors as large as ~F ULPs have been seen."
                 diff))
    (pnl (format nil "y3n errors as large as ~F ULPs."
                 adiffx3n))
    (pnl (format nil "yyyn error as large as ~F ULPs."
                 adiffxxxn))
    ))
