;;; $Id: getsurf.scm,v 1.1 2005/01/16 20:41:03 mstorti Exp $
(load "./femref.scm")
#!
(getsurf "./cube.con.tmp" "./cube.nod.tmp" "./cube.state.tmp" 
	 "cube.surf-con.tmp" "cube.grad-u.tmp" 1)
!#

#!
(define v (make-dvdbl))
(define N 100)
(dvdbl-resize! v N)
(define j 0)
(while 
 (< j N)
 (dvdbl-set! v j j)
 (set! j (+ j 1)))
!#

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
	    surf-mass node-mass u us)
