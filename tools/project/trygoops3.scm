(use-modules (oop goops))
(use-modules (ice-9 format))

;;; The generic function print in base class
;;; uses the functions of the derived A and B.

(define start 0)
(define (l-init)
;  (format #t "start ~A\n" start)
  (let loop ((j 0)
	     (q '()))
    (cond ((= j 5) 
	   (set! start (+ start j))
	   (format #t "new l-part ~A\n" q)
	   (reverse q))
	  (else 
	   (loop (+ 1 j) (cons (+ j start) q))))))

(define-class <A>()
  (l #:accessor l-part))

(define-method (initialize (a <A>) . args)
  (set! (l-part a) (l-init)))

(define a1 (make <A>))
(define a2 (make <A>))

(format #t "a1.l ~A, a2.l ~A\n" (l-part a1) (l-part a2))

(define a3 (make (class-of a1)))
(format #t "a3 ~A\n" a3)

