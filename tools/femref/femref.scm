;(load-extension "./libfemref" "init_femref")
;(define (my-double . args) (format #t "args ~A\n" args) (* 2 x))
;(getsurf my-double)

(define (my-double x) (* 2 x))

(format #t "2*2 = ~A\n" (apply my-double 2))
