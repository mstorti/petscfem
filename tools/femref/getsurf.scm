;;; $Id: getsurf.scm,v 1.17 2005/01/18 20:38:17 mstorti Exp $
(load "./femref.scm")
(use-modules (dvector))
(load "./utils.scm")

(define ctx (make-get-surf-ctx))
(define icone (make-dvint))
(define surf-con (make-dvint))
(define surf-nodes (make-dvint))
(define x (make-dvdbl))
(define surf-mass (make-dvdbl))
(define node-mass (make-dvdbl))
(define ndof 3)
(define ndim 3)

(format #t "icone: read ~A ints\n" (dvint-cat! icone "cube.con.tmp"))
(dvdbl-reshape! x 0 ndim)
(format #t "xnod: read ~A dbls\n" (dvdbl-cat! x "cube.nod.tmp"))

(getsurf ctx icone surf-con surf-nodes 1 0)
(define nfaces (car (dvint-shape surf-con)))

(comp-matrices ctx surf-con
	       surf-nodes x surf-mass node-mass)
(define nnod (car (dvdbl-shape x)))
(define nsurf-nodes (dvint-size surf-nodes))

(define un (make-dvdbl))
(dvdbl-resize! un nsurf-nodes ndof)

(define ue (make-dvdbl))
(dvdbl-resize! ue nfaces ndof)
(define uee (make-dvdbl))
(dvdbl-clone! uee ue)

(nod->elem-proj ctx surf-con un uee)
;(ddump "uee" uee)

(elem->nod-proj ctx surf-con surf-mass node-mass uee un)
;(ddump "un" un)
