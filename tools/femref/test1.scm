;;; $Id: test1.scm,v 1.1 2005/01/18 15:23:29 mstorti Exp $
(load "./femref.scm")
(use-modules (dvector))
(load "./utils.scm")

(define x (make-dvdbl))
(dvdbl-resize! x 10 3)
(dvdbl-apply! x (lambda (x) (- (* 2 (random 10)))))

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

(dvdbl-set! x (lambda (indx) (+ (* (car indx) 1000) (cadr indx))))
(ddump "x = 1000*j+k" x)
