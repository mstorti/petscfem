;;; $Id: test3.scm,v 1.1 2005/02/07 13:13:27 mstorti Exp $
(load "./femref.scm")
(use-modules (dvector))
(load "./utils.scm")

(define x (make-dvdbl))
(dvdbl-resize! x 10 3)
(dvdbl-apply! x (lambda (x) (- (* 2 (random 10)))))

(format #t "x ~A\n" x)
