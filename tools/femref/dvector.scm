;;; $Id: dvector.scm,v 1.4 2005/01/18 15:23:29 mstorti Exp $
(define-module (dvector))

(load-extension "./libfemref" "dvint_init")
(load-extension "./libfemref" "dvdbl_init")

(define-public (dvdbl-resize! v . shape)
  (let loop ((size 1)
	     (q shape))
    (cond ((null? q) 
	   (format #t "size ~A, shape ~A\n" size shape)
	   (dvdbl-resize-w! v size)
	   (apply dvdbl-reshape! v shape))
	  (#t (loop (* size (car q)) (cdr q))))))
	   
(define (dvdbl-set-with-filler v filler)
  (let ((shape (dvdbl-shape v)))
    (let loop ((q shape)
	       (indx '()))
      (cond ((null? q) 
	     (dvdbl-set! v (reverse indx) 
			 (filler (reverse indx))))
	    (#t (let ((n (car q)))
		  (do ((j 0 (+ j 1))) ((= j n))
		    (loop (cdr q) (cons j indx)))))))))

(define-public (dvdbl-set! v . args)
  (cond ((= (length args) 2) (apply dvdbl-set-w2 v args))
	(#t (let ((arg (car args)))
	      (cond ((procedure? arg) (dvdbl-set-with-filler v arg))
		    (#t (apply dvdbl-set-w1 v args)))))))

(define-public (dvdbl-apply! v fun)
  (let ((n (dvdbl-size v)))
    (do ((j 0 (+ j 1))) ((= j n)) 
      (let ((w (dvdbl-ref v j)))
	(dvdbl-set! v j (fun w))))))

(define-public (dvdbl-add! v alpha)
  (dvdbl-apply! v (lambda (y) (+ alpha y))))

(define-public (dvdbl-rand! v)
  (dvdbl-apply! v (lambda (y) (random:uniform))))

(define-public (dvdbl-assoc v fun init)
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

(define-public (dvdbl-max v . args)
  (cond ((null? args) (dvdbl-max-aux v >))
	(#t (dvdbl-max-aux v (car args)))))
		
(define (dvdbl-min-aux v comp)
  (dvdbl-max v (lambda (x y) (not (comp x y)))))

(define-public (dvdbl-min v . args)
  (cond ((null? args) (dvdbl-max-aux v <))
	(#t (dvdbl-min-aux v (car args)))))

(export make-dvint dvint-clone! dvint-push! dvint-size
	dvint-reshape! dvint-shape dvint-ref
	dvint-read! dvint-cat! dvint-dump)

(export make-dvdbl dvdbl-clone! dvdbl-push! dvdbl-size 
	dvdbl-reshape! dvdbl-shape dvdbl-ref
	dvdbl-read! dvdbl-cat! dvdbl-dump)

(export dvdbl-scale!)
