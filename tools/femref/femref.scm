;;; $Id: femref.scm,v 1.19 2005/01/17 03:45:45 mstorti Exp $
(use-modules (ice-9 optargs))
(load-extension "./libfemref" "init_femref")
(load "./while2.scm")

(define* (fem-smooth ctx surf-con elem-mass node-mass 
		     u us #:key (niter 5) (verbose #f))
  (fem-smooth-w ctx surf-con elem-mass node-mass 
		u us niter verbose))
