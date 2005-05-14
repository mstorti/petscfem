(set! %load-path (cons "../femref" %load-path))

(use-modules (oop goops))
(use-modules (ice-9 format))
(load-from-path "utils.scm")
(load-from-path "while2")
(use-modules (dvector2))

(define v (make <dvdbl>))
(dv-resize! v 2 3)

(dv-dump v)
;;(dvdbl-dump (vec v))
