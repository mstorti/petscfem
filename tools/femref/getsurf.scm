;;; $Id: getsurf.scm,v 1.16 2005/01/18 15:23:29 mstorti Exp $
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

(ddump "x: " x)
(dvdbl-scale! x 4.)
(ddump "x *= 4" x)

(define (dvdbl-scaleb! v alpha)
  (dvdbl-apply! v (lambda (y) (* alpha y))))

(dvdbl-scaleb! x 0.5)
(ddump "x /= 2" x)
(dvdbl-add! x 100)
(ddump "x += 100" x)
(format #t "max x ~A\n"  (dvdbl-max x))
(format #t "min x ~A\n"  (dvdbl-min x))
(dvdbl-rand! x )
(dvdbl-apply! x (lambda (x) (- (* 2 x) 1)))
(ddump "x = rand" x)
(define M 1000)
(dvdbl-apply! x (lambda (x) (- (* 2 (random M)))))
(ddump "x = (rand 1000)" x)
(format #t "max x ~A\n"  (dvdbl-max x))
(format #t "max |x| ~A\n"  
	(dvdbl-assoc x (lambda (x y) 
			 (cond ((> ($abs x) ($abs y)) x)
			       (#t y))) 0))
(quit)

(getsurf ctx icone surf-con surf-nodes 1 0)
(define nfaces (car (dvint-shape surf-con)))

(comp-matrices ctx surf-con
	       surf-nodes x surf-mass node-mass)
(define nnod (car (dvdbl-shape x)))
(define nsurf-nodes (dvint-size surf-nodes))

(define un (make-dvdbl))
(dvdbl-resize! un nsurf-nodes ndof)

(do ((j 0 (+ j 1))) ((= j nsurf-nodes)) 
  (let ((node (dvint-ref surf-nodes j)))
    (do ((k 0 (+ k 1))) ((= k ndim))
      (dvdbl-set! un (list j k) (dvdbl-ref x (list node k))))))
(ddump "un" un)

(define ue (make-dvdbl))
(dvdbl-resize! ue nfaces ndof)
(define uee (make-dvdbl))
(dvdbl-clone! uee ue)

(nod->elem-proj ctx surf-con un uee)
(ddump "uee" uee)

(elem->nod-proj ctx surf-con surf-mass node-mass uee un)
(ddump "un" un)
