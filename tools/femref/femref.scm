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

(define icone (make-dvint))
(format #t "read ~A ints\n" (dvint-cat! icone "cube.con.tmp"))
(getsurf icone 5)
