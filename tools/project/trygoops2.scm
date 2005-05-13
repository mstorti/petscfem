(use-modules (oop goops))
(use-modules (ice-9 format))

;;; The generic function print in base class
;;; uses the functions of the derived A and B.

(define-class super ())
(define-method (print (s super) v) 
  (let loop ((q v))
    (cond ((not (null? q)) 
	   (print-elem s (car q))
	   (loop (cdr q))))))

(define-class A (super))
(define-method (meth (a A)) (format #t "in meth A\n"))
(define-method (print-elem (a A) x) (format #t "A element ~s\n" x))

(define-class B (super))
(define-method (meth (b B)) (format #t "in meth B\n"))
(define-method (print-elem (b B) x) (format #t "B element ~s\n" x))

(define a (make A))
(define b (make B))

(print a '(1 2 3 4 5))
(print b '(1 2 3 4 5))
