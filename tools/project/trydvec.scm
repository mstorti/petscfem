(set! %load-path (cons "../femref" %load-path))

(use-modules (oop goops))
(use-modules (ice-9 format))
(load-from-path "utils.scm")
(load-from-path "while2")
(use-modules (dvector))

(define v (make <dvdbl>))
(define w (make <dvdbl>))

(dv-resize! v 5 4)
;(define (dv-dump-v) (dv-dump v) (newline))
(dv-dump v "v size 5x4")

(dv-set! v 
  (lambda (indx)
    (let loop ((q indx)
	       (x 0))
      (cond ((null? q) x)
	    (else (loop (cdr q) (+ (* x 100) (car q))))))))
(dv-dump v "v filled 100*j+k")

(define (m-set v M)
  (dv-set! v
	   (lambda (indx)
	     (let loop ((q indx)
			(x 0))
	       (cond ((null? q) x)
		     (else (loop (cdr q) (+ (* x M) (car q)))))))))
(m-set v 1000)
(dv-dump v "v filled 1000*j+k")

(dv-set! v 7)
(dv-dump v "v <- 7")

(dv-set! v (lambda (indx) (+ (cadr indx) 0.5)))
(dv-dump v "v <- j+0.5")

(dv-set! v (lambda (indx) (car indx)))
(dv-dump v "v <- k")

(define w (make <dvdbl>))

(dv-clone! w v)
(m-set w 1000)
(dv-slice-range! v w '(0 0 5 2))
(dv-dump v "v <- slice ~A\n" '(0 0 5 2))

(dv-slice! v w '(0 #f #f 2) '(1 #f #f 2))
(dv-dump v "v <- slice ~A\n" '((0 #f #f 2) (1 #f #f 2)))

(define x (make <dvint>))
(dv-resize! x 7 6)

(dv-set! x (lambda (indx) (inexact->exact (* 23.24253 (car indx)))))
(dv-dump x "x integer vector")
