;;; $Id: getsurf.scm,v 1.7 2005/01/17 15:43:26 mstorti Exp $
(load "./dvector.scm")
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
(define nfaces (car (dvint-shape surf-con)))

(comp-matrices ctx surf-con
	       surf-nodes x surf-mass node-mass)
(define nnod (car (dvdbl-shape x)))
(define nsurf-nodes (dvint-size surf-nodes))

(define ndof 3)
(define ue (make-dvdbl))
(dvdbl-resize! ue nfaces ndof)
(dvdbl-set! ue 5)
(define un (make-dvdbl))
(dvdbl-resize! un nsurf-nodes ndof)

(elem->nod-proj ctx surf-con surf-mass node-mass ue un)
(format #t "\n\nun:\n")
(dvdbl-dump un)

; (fem-smooth ctx surf-con surf-nodes
; 	    surf-mass node-mass u us #:verbose #f)
