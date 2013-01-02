;; ran.lisp -- random number generator from Numerical Recipes


(defmacro dfloat (x)
  `(float ,x 1d0))

(defvar ran1 

(let 
 ((iff 0)
  (ix1 0)
  (ix2 0)
  (ix3 0)
  (r (make-array 97 :element-type 'double-float :initial-element 0d0)))
 (declare (type fixnum iff ix1 ix2 ix3))
 (declare (type (simple-array double-float (*)) r))

#'(lambda (idum)
   (declare (type fixnum idum))

   (prog ((m1 0) (ia1 0) (ic1 0) (m2 0) (ia2 0) (ic2 0) (rm1 0d0) (rm2 0d0)
          (m3 0) (ia3 0) (ic3 0) (j 0) (ran1 0d0))
    (declare (type fixnum m1 ia1 ic1 m2 ia2 ic2 m3 ia3 ic3 j))
    (declare (type double-float rm1 rm2 ran1))

  (setq m1 259200
        ia1 7141
        ic1 54773
        rm1 3.8580247d-6) 
  (setq m2 134456
        ia2 8121
        ic2 28411
        rm2 7.4373773d-6) 
  (setq m3 243000
        ia3 4561
        ic3 51349) 

  (when 
   (or (< idum 0) (= iff 0)) 
   (setf iff 1)
   (setf ix1 (mod (- ic1 idum) m1))
   (setf ix1 (mod (+ (* ia1 ix1) ic1) m1)) 
   (setf ix2 (mod ix1 m2))
   (setf ix1 (mod (+ (* ia1 ix1) ic1) m1)) 
   (setf ix3 (mod ix1 m3))

   (do ((j 1 (1+ j)))
       ((> j 97) t)
       (declare (type fixnum j))
     (setf ix1 (mod (+ (* ia1 ix1) ic1) m1))
     (setf ix2 (mod (+ (* ia2 ix2) ic2) m2))
     (setf (aref r (1- j)) (* (+ (dfloat ix1) (* (dfloat ix2) rm2)) rm1)))
   (setf idum 1)) 

  (setf ix1 (mod (+ (* ia1 ix1) ic1) m1)) 
  (setf ix2 (mod (+ (* ia2 ix2) ic2) m2)) 
  (setf ix3 (mod (+ (* ia3 ix3) ic3) m3)) 
  (setf j (+ 1 (floor (/ (* 97 ix3) m3)))) 
  (if (or (> j 97) (< j 1)) (error " error in ran1 ")) 
  (setf ran1 (aref r (1- j))) 
  (setf (aref r (1- j)) (* (+ (dfloat ix1) (* (dfloat ix2) rm2)) rm1)) 
   
  (return (the double-float ran1))))))
