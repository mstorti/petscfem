(load-extension "./libfemref" "init_femref")
(define (my-double x) (* 2 x))
; (define (my-double . x) (format #t "#args: ~A\n" (length x)))
(getsurf my-double)

#!
(define (my-double x) (* 2 x))
(format #t "2*2 = ~A\n" (apply my-double 2 '()))
!#
