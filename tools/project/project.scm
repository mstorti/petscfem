;;; $Id: project.scm,v 1.1 2005/03/03 02:21:00 mstorti Exp $
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
  (dvint-set! indx j (+ j 1)))
(idump "indx" indx)

(dvdbl-slice-indx! w x indx 0)
