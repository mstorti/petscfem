;;; $Id: dvector.scm,v 1.2 2005/01/18 02:29:47 mstorti Exp $

(load-extension "./libfemref" "dvint_init")
(load-extension "./libfemref" "dvdbl_init")

(define (dvdbl-resize! v . shape)
  (let loop ((size 1)
	     (q shape))
    (cond ((null? q) 
	   (format #t "size ~A, shape ~A\n" size shape)
	   (dvdbl-resize-w! v size)
	   (apply dvdbl-reshape! v shape))
	  (#t (loop (* size (car q)) (cdr q))))))
	   
(define (dvdbl-set! v . args)
  (cond ((= (length args) 2) (apply dvdbl-set-w2 v args))
	(#t (apply dvdbl-set-w1 v args))))

(define (dvdbl-apply! v fun)
  (let ((n (dvdbl-size v)))
    (do ((j 0 (+ j 1))) ((= j n)) 
      (let ((w (dvdbl-ref v j)))
	(dvdbl-set! v j (fun w))))))

(define (dvdbl-scaleb! v alpha)
  (dvdbl-apply! v (lambda (y) (* alpha y))))

(define (dvdbl-add! v alpha)
  (dvdbl-apply! v (lambda (y) (+ alpha y))))

(define (dvdbl-rand! v)
  (dvdbl-apply! v (lambda (y) (random:uniform))))

(define (dvdbl-assoc v fun init)
  (let ((n (dvdbl-size v)))
    (cond ((= n 0) init)
	  (#t (let ((result (dvdbl-ref v 0)))
		(do ((j 1 (+ j 1))) ((= j n))
		  (set! result (fun (dvdbl-ref v j) result)))
		result)))))

(define (dvdbl-max-aux v comp)
  (dvdbl-assoc v 
	       (lambda (x y) 
		 (cond ((comp x y) x) (#t y)))
	       0))

(define (dvdbl-max v . args)
  (cond ((null? args) (dvdbl-max-aux v >))
	(#t (dvdbl-max-aux v (car args)))))
		
(define (dvdbl-min-aux v comp)
  (dvdbl-max v (lambda (x y) (not (comp x y)))))

(define (dvdbl-min v . args)
  (cond ((null? args) (dvdbl-max-aux v <))
	(#t (dvdbl-min-aux v (car args)))))
