;;; $Id: femref.scm,v 1.17 2005/01/16 21:14:36 mstorti Exp $
(use-modules (ice-9 optargs))
(load-extension "./libfemref" "init_femref")
(load-extension "./libfemref" "dvint_init")
(load-extension "./libfemref" "dvdbl_init")
(load "./while2.scm")

(define* (fem-smooth ctx surf-con surf-nodes 
		     surf-mass node-mass u us #:key (niter 5) (verbose #f))
  (fem-smooth-w ctx surf-con surf-nodes 
		surf-mass node-mass u us niter verbose))
