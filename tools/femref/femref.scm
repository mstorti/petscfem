;;; $Id: femref.scm,v 1.20 2005/03/03 02:20:56 mstorti Exp $
(use-modules (ice-9 optargs))
(load-extension "../femref/libfemref" "init_femref")
(load-from-path "while2.scm")

(define* (fem-smooth ctx surf-con elem-mass node-mass 
		     u us #:key (niter 5) (verbose #f))
  (fem-smooth-w ctx surf-con elem-mass node-mass 
		u us niter verbose))
