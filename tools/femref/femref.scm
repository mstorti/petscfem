(load-extension "./libfemref" "init_femref")
(load-extension "./libfemref" "dvint_init")
(load-extension "./libfemref" "dvdbl_init")
(load "./while2.scm")

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

(define get-surf-ctx (make-get-surf-ctx))
(define icone (make-dvint))
(define surf-con (make-dvint))
(define surf-nodes (make-dvint))
(format #t "read ~A ints\n" (dvint-cat! icone "cube.con.tmp"))

(getsurf get-surf-ctx icone surf-con surf-nodes 1 0)
(format #t "surf-con:\n")
(dvint-dump surf-con)

(format #t "surf-nodes:\n")
(dvint-dump surf-nodes)
