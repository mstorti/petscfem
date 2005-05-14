;;; $Id: dvector2.scm,v 1.1 2005/05/14 16:01:17 mstorti Exp $
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

(define-macro (dv-method type fun)
  `(define-method (,(dv-fun fun) (v ,(dv-class type)) . rest)
     (apply ,(dvtype-fun type fun) (vec v) rest)))

(define-class <dvdbl> (<dvector>)
  (v #:init-value (make-dvdbl)
     #:accessor vec))

#!
(define-method (dv-resize-w! (v <dvdbl>) size)
  (dvdbl-resize-w! (vec v) size))
!#

(dv-method dvdbl resize-w!)

(define-method (dv-reshape! (v <dvdbl>) . shape)
  (apply dvdbl-reshape! (vec v) shape))

(define-method (dv-dump (v <dvdbl>) . rest)
  (apply dvdbl-dump (vec v) rest))

#!
(defmacro (define-dv-method class type fun)
  `(define-method (dv-dump (v ,class) . rest)
     (apply dvdbl-dump (vec v) rest)))

  scm_c_define_gsubr(DVTYPE "-clone!", 2, 0, 0, scm_fun(DVECTOR_CLONE_FUN));
  scm_c_define_gsubr(DVTYPE "-push!", 2, 0, 0, scm_fun(DVECTOR_PUSH_FUN));
  scm_c_define_gsubr(DVTYPE "-size", 1, 0, 0, scm_fun(DVECTOR_SIZE_FUN));
  scm_c_define_gsubr(DVTYPE "-resize-w!", 2, 0, 0, scm_fun(DVECTOR_RESIZE_FUN));
  scm_c_define_gsubr(DVTYPE "-reshape!", 1, 0, 1, scm_fun(DVECTOR_RESHAPE_FUN));
  scm_c_define_gsubr(DVTYPE "-shape", 1, 0, 0, scm_fun(DVECTOR_SHAPE_FUN));
  scm_c_define_gsubr(DVTYPE "-set-w1", 2, 0, 0, scm_fun(DVECTOR_SET_W1_FUN));
  scm_c_define_gsubr(DVTYPE "-set-w2", 3, 0, 0, scm_fun(DVECTOR_SET_W2_FUN));
  scm_c_define_gsubr(DVTYPE "-ref", 2, 0, 0, scm_fun(DVECTOR_REF_FUN));
  scm_c_define_gsubr(DVTYPE "-read!", 2, 0, 0, scm_fun(DVECTOR_READ_FUN));
  scm_c_define_gsubr(DVTYPE "-cat!", 2, 0, 0, scm_fun(DVECTOR_CAT_FUN));
  scm_c_define_gsubr(DVTYPE "-dump", 1, 2, 0, scm_fun(DVECTOR_DUMP_FUN));
!#

(export <dvector> <dvdbl> dvdbl-dump 
	dv-dump dv-reshape! dv-resize! vec)
