(load-extension "./libfemref" "init_femref")

#!
(define (my-double x) (* 2 x))
(getsurf my-double)
!#
  
(getsurf "./cube.con.tmp" "./cube.nod.tmp" "./cube.state.tmp" 
	 "cube.surf-con.tmp" "cube.grad-u.tmp" 1)
