;;; $Id: dvector2.scm,v 1.2 2005/05/14 16:36:33 mstorti Exp $
(define-module (dvector2))
(use-modules (oop goops))

(load-extension 
 (search-path %load-path "libfemref.so") 
 "dvint_init")
(load-extension 
 (search-path %load-path "libfemref.so") 
 "dvdbl_init")

(define-class <dvector>())

(define-method (dv-resize! (v <dvector>) . shape)
  (let loop ((size 1)
	     (q shape))
    (cond ((null? q) 
	   (dv-resize-w! v size)
	   (apply dv-reshape! v shape))
	  (#t (loop (* size (car q)) (cdr q))))))

(define-macro (dv-class type)
  `(string->symbol (string-append "<" (symbol->string ,type) ">")))

(define-macro (dv-fun fun)
  `(string->symbol (string-append 
		    "dv-"
		    (symbol->string ,fun))))

(define-macro (dvtype-fun type fun)
  `(string->symbol (string-append 
		    (symbol->string ,type) "-"
		    (symbol->string ,fun))))

(define-macro (dv-method1 type fun)
  `(define-method (,(dv-fun fun) (v ,(dv-class type)) . rest)
     (apply ,(dvtype-fun type fun) (vec v) rest)))

(define-macro (dv-method fun)
  `(dv-method1 dvdbl ,fun)
  `(dv-method1 dvint ,fun))

; (define-macro (dv-method fun)
;   `(dv-method1 dvdbl ,fun)
;   `(dv-method1 dvint ,fun)
;   `(export ,(dv-fun fun)))

(define-class <dvdbl> (<dvector>)
  (v #:init-value (make-dvdbl)
     #:accessor vec))

(dv-method resize-w!)
(dv-method clone!)
(dv-method push!)
(dv-method size)
(dv-method reshape!)
(dv-method shape)
(dv-method set-w1)
(dv-method set-w2)
(dv-method ref)
(dv-method read!)
(dv-method cat!)
(dv-method dump)

(export dv-resize-w!)
(export dv-clone!)
(export dv-push!)
(export dv-size)
(export dv-reshape!)
(export dv-shape)
(export dv-set-w1)
(export dv-set-w2)
(export dv-ref)
(export dv-read!)
(export dv-cat!)
(export dv-dump)

(export <dvector> <dvdbl> vec)
