(set! %load-path (cons "../femref" %load-path))

(use-modules (oop goops))
(use-modules (ice-9 format))
(load-from-path "utils.scm")
(load-from-path "while2")
(use-modules (dvector2))

(define v (make <dvdbl>))
(dv-resize! v 2 3)

(dv-dump v)

(dv-set-with-filler! v 
;  (lambda (j k) (+ (* 1000 j) k))
  (lambda (indx)
    (let loop ((q indx)
	       (x 0))
      (cond ((null? q) x)
	    (else (loop (cdr q) (+ (* x 100) (car q))))))))
(dv-dump v)
