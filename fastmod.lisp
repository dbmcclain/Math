;;; Fast  (mod (expt number exponent) modulus)
(defun expt-mod (number exponent modulus)
  (loop with result = 1
        for exp = exponent then (floor exp 2)
        for sqr = number then (mod (* sqr sqr) modulus)
        until (zerop exp)
        when (oddp exp) do
        (setf result (mod (* result sqr) modulus))
        finally (return result)))

;;Nothing special so far...

;;Here's some more...

;;; if d = gcd(a,b) then there is a x and b y so that d = ax + by
(defun extended-gcd (a b)
  (loop for x = 1 then x2
        and y = 0 then y2
        and d = a then d2
        and x2 = 0 then (- x (* k x2))
        and y2 = 1 then (- y (* k y2))
        and d2 = b then (- d (* k d2))
        until (zerop d2)
        for k = (floor d d2)
        finally (return (values d x y))))

;;;  a^-1 mod n , provided that it exists
(defun multiplicative-inverse-mod (a n)
  (multiple-value-bind (d x y)
               (extended-gcd a n)
               (if (> d 1)
                   (error "INVERSE-MOD: no inverse~
                     => modulus is not prime.")
                 (loop while (minusp x) do
                   (incf x n)
                   finally (return x)))))

;I've used it for implementing a miller-rabin pseudo-prime test and an 
;implementation of the RSA public key cipher.
;If I find the time I will place this code (and the rest like miller-rabin 
;and RSA) on my Webpage http://www.dataheaven.de.
;Meanwhile permission is hereby granted to use the above code if you include 
;the following copyright statement and disclaimer. (I know it's not very 
;much code but it's part of some more code that I worked some time to get it 
;working as it does)

;;;;                                                                            ;
;;;; (c) 2001 by Jochen Schmidt.
;;;;
;;;; Redistribution and use in source and binary forms, with or without
;;;; modification, are permitted provided that the following conditions
;;;; are met:
;;;; 1. Redistributions of source code must retain the above copyright
;;;;    notice, this list of conditions and the following disclaimer.
;;;; 2. Redistributions in binary form must reproduce the above copyright
;;;;    notice, this list of conditions and the following disclaimer in the
;;;;    documentation and/or other materials provided with the distribution.
;;;;
;;;; THIS SOFTWARE IS PROVIDED "AS IS" AND THERE ARE NEITHER 
;;;; EXPRESSED NOR IMPLIED WARRANTIES -  THIS INCLUDES, BUT 
;;;; IS NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
;;;; AND FITNESS FOR A PARTICULAR PURPOSE.IN NO WAY ARE THE
;;;; AUTHORS LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;;; SPECIAL, EXEMPLARY OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
;;;; NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES ;
;;;; LOSS OF USE, DATA, OR PROFITS ; OR BUSINESS INTERRUPTION)
;;;; 
;;;; For further details contact the authors of this software.
;;;;
;;;;  Jochen Schmidt        
;;;;  Zuckmantelstr. 11     
;;;;  91616 Neusitz         
;;;;  GERMANY               
;;;;
;;;; Nürnberg, 16.Apr.2001 Jochen Schmidt
;;;;

