(use-modules (oop goops))
(use-modules (ice-9 format))

;;; Example of use of classes

(define-class A ())
(define-method (meth (a A)) (format #t "in meth A\n"))

(define-class B ())
(define-method (meth (b B)) (format #t "in meth B\n"))

(define a (make A))
(define b (make B))

(meth a)
(meth b)
