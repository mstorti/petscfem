;;; $Id: trykeyw.scm,v 1.1 2005/05/17 16:29:39 mstorti Exp $
(use-modules (ice-9 format))
(use-modules (ice-9 optargs))

(define* (f x y #:key a b . rest)
	   (format #t "x ~A, y ~A, a ~A, b ~A, rest ~A\n"
		   x y a b rest))

(f 1 2 3)
(f 1 2 3 #:a 4)
(f 1 2 #:a 4 3)

(define* (g #:key a b)
	   (format #t "a ~A, b ~A\n" a b))

(g #:a 1 #:b 2)
(g #:b 1 #:a 2)
