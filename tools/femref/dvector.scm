;;; $Id: dvector.scm,v 1.1 2005/01/17 03:45:45 mstorti Exp $

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
