;;; $Id: getsurf.scm,v 1.2 2005/01/16 21:14:36 mstorti Exp $
(load "./femref.scm")

(define ctx (make-get-surf-ctx))
(define icone (make-dvint))
(define surf-con (make-dvint))
(define surf-nodes (make-dvint))
(define x (make-dvdbl))
(define surf-mass (make-dvdbl))
(define node-mass (make-dvdbl))

(format #t "icone: read ~A ints\n" (dvint-cat! icone "cube.con.tmp"))
(format #t "xnod: read ~A dbls\n" (dvdbl-cat! x "cube.nod.tmp"))

(getsurf ctx icone surf-con surf-nodes 1 0)
(format #t "surf-con:\n")
(dvint-dump surf-con)

(format #t "surf-nodes:\n")
(dvint-dump surf-nodes)

(comp-matrices ctx surf-con
	       surf-nodes x surf-mass node-mass)

(format #t "\n\nsurf-mass:\n")
(dvdbl-dump surf-mass)

(format #t "\n\nnode-mass:\n")
(dvdbl-dump node-mass)

(define ndof 4)
(define u (make-dvdbl))
(define us (make-dvdbl))
(dvdbl-cat! u "cube.state.tmp")

(fem-smooth ctx surf-con surf-nodes
	    surf-mass node-mass u us #:verbose #f)
