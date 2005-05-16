;;; $Id: trydvec.scm,v 1.20 2005/05/16 03:25:08 mstorti Exp $
(set! %load-path (cons "../femref" %load-path))

(use-modules (oop goops))
(use-modules (ice-9 format))
(load-from-path "utils.scm")
(load-from-path "while2")
(use-modules (dvector))

(define v (make <dvdbl>))
(define w (make <dvdbl>))
(define (m-set v M)
  (dv-set! v
	   (lambda (indx)
	     (let loop ((q indx)
			(x 0))
	       (cond ((null? q) x)
		     (else (loop (cdr q) (+ (* x M) (car q)))))))))
(define ivec (make <dvint>))
(define x (make <dvint>))
(define w2 (make <dvdbl>))
(define w3 (make <dvint>))
(define w4 (make <dvint>))
(define (my-sum . args) 
  (let ((n (length args)))
    (cond ((= n 0) 0)
	  ((= n 1) (apply + args)))))

(dv-resize! v 5 4)
(dv-dump v "qqqq" 1 "v size 5x4")

#!
(dv-set! v 
	 (lambda (indx)
	   (let loop ((q indx)
		      (x 0))
	     (cond ((null? q) x)
		   (else (loop (cdr q) (+ (* x 100) (car q))))))))
(dv-dump v "v filled 100*j+k")

(m-set v 1000)
(dv-dump v "v filled 1000*j+k")

(dv-set! v 7)
(dv-dump v "v <- 7")

(dv-set! v (lambda (indx) (+ (cadr indx) 0.5)))
(dv-dump v "v <- j+0.5")

(dv-set! v (lambda (indx) (car indx)))
(dv-dump v "v <- k")

(dv-clone! w v)
(m-set w 1000)
(dv-slice-range! v w '(0 5 2) 0)
(dv-dump v "v <- slice ~A\n" '(0 0 5 2))

(dv-slice! v w '(#f #f 2)  0 '(#f #f 2) 1)
(dv-dump v "v <- slice ~A\n" '((0 #f #f 2) (1 #f #f 2)))

(dv-resize! ivec 3)
(dv-set! ivec (lambda(q) (car q)))
(dv-dump ivec "ivec: ")

(dv-slice! v w ivec 1)
(dv-dump v "w(1:3,:) ")

(dv-resize! x 7 6)

(dv-set! x (lambda (indx) (inexact->exact (* 23.4 (car indx)))))
(dv-dump x "x integer vector")

(dv-resize! w2 5 5)
(m-set w2 1000)
(dv-dump w2 "w2 filled with 1000*j+k")

(dv-apply! w2 (lambda(x) (sqrt (+ x 1))))
(dv-dump w2 "w2 <- sqrt(w2+1)")
(set! w2 #f)

(dv-resize! w3 5 5)
(m-set w3 10)
(dv-dump w3 "w3 filled with 10*j+k")

(dv-clone! w4 w3)

(dv-apply! w3 (lambda(x) (modulo x 10)))
(dv-dump w3 "w3 <- (modulo x 10)")
(set! w3 #f)

(dv-apply! w4 (lambda(x) (- (* 2 (random 10)) 10)))
(dv-dump w4 "w4 <- random [-5 5]")
(format #t "sum w4: ~A\n" (dv-reduce w4 (lambda(x y) (+ x y)) 0))
(format #t "sum(|w4|): ~A\n" (dv-reduce w4 (lambda(x y) (+ (abs x) (abs y))) 0))

(set! w4 (make <dvdbl>))
(format #t "sum w4 (empty): ~A\n" (dv-reduce w4 my-sum))

(dv-resize! w4 1)
(dv-set! w4 23)
(format #t "sum w4 [23]: ~A\n" (dv-reduce w4 my-sum))

(format #t "version ~A\n" dv-version)

(dv-resize! w4 7 7)
(format #t "antes de dv-rand!\n")
(dv-rand! w4)
(format #t "despues de dv-rand!\n")

(dv-dump w4 "w4 dbl rand (rand!)")

(set! w4 (make <dvint>))
(dv-resize! w4 7 7)
(dv-rand! w4 100)
(dv-dump w4 "w4 int rand (rand!)")

; #!
; (define w5 (make <dvdbl>))
; (dv-resize! w5 6 6)
; (dv-rand! w5)
; (dv-dump w5 "w5 rand" 1)

; (set! w4 (make <dvint>))
; (dv-resize! w4 6 6)
; (dv-rand! w4 10)
; (dv-dump w4 "w4 rand")
; !#
!#

