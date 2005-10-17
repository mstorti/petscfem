;;; $Id: project.scm,v 1.2 2005/10/17 02:51:37 mstorti Exp $
(set! %load-path (cons "../femref" %load-path))

;(use-modules (femref))
(load-from-path "utils.scm")
(load-from-path "while2")
(use-modules (dvector))

(load-extension "./libproject" "init_fem_interp")

(define fem-interp (make-fem-interp))
(format #t "~A" fem-interp)
(set! fem-interp #f)

(define v (make-dvint))
(do ((j 0 (+ j 1))) ((= j 4))
  (dv-set! indx j (+ j 1)))
(idump "indx" indx)

(dv-slice-indx! w x indx 0)
