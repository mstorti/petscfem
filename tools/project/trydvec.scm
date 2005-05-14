(set! %load-path (cons "../femref" %load-path))

(use-modules (oop goops))
(use-modules (ice-9 format))
(load-from-path "utils.scm")
(load-from-path "while2")
(use-modules (dvector))

(define v (make <dvdbl>))
(dv-resize! v 5 4)
(define (dv-dump-v) (dv-dump v) (newline))
(dv-dump-v)

(dv-set! v 
  (lambda (indx)
    (let loop ((q indx)
	       (x 0))
      (cond ((null? q) x)
	    (else (loop (cdr q) (+ (* x 100) (car q))))))))
(dv-dump-v)

(dv-set! v
  (lambda (indx)
    (let loop ((q indx)
	       (x 0))
      (cond ((null? q) x)
	    (else (loop (cdr q) (+ (* x 1000) (car q))))))))
(dv-dump-v)

(dv-set! v 7)
(dv-dump-v)

(dv-set! v (lambda (indx) (car indx)))
(dv-dump-v)

(dv-set! v (lambda (indx) (cadr indx)))

(define w (make <dvdbl>))
(format #t "v size ~A\n" (dv-size v))
(format #t "w size ~A\n" (dv-size w))

(dv-clone! w v)

(dv-slice-range! w v '(1 1 5 2))
(dv-dump-v)
