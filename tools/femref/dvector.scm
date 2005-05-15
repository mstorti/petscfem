;;; $Id: dvector.scm,v 1.9 2005/05/15 13:28:12 mstorti Exp $
(define-module (dvector))
(use-modules (oop goops))

(load-extension 
 (search-path %load-path "libfemref.so") 
 "dvint_init")
(load-extension 
 (search-path %load-path "libfemref.so") 
 "dvdbl_init")

;;; We first define the base class <dvector>
;;; and then add some functions that operate
;;; the same way on <dvdbl> and <dvint>
(define-class <dvector>())

;;; Resizes and reshapes a vector to new shape `shape'
;;; For example: (dv-resize! v 2 3 4)
(define-method (dv-resize! (v <dvector>) . shape)
  (let loop ((size 1)
	     (q shape))
    (cond ((null? q) 
	   (dv-resize-w! v size)
	   (apply dv-reshape! v shape))
	  (#t (loop (* size (car q)) (cdr q))))))

;;; Fills V with filler function `filler'.
;;; It's equivalent to set each element `(j k l ...)'
;;; with (filler j k l ...)
;;; Examples: 
;;; (dv-set-with-filler! v (lambda (indx) 7))
;;; sets all elements to 7.
;;; (dv-set-with-filler! v (lambda (indx) (car indx)))
;;; sets all elements to the number of row.
(define-public (dv-set-with-filler! v filler)
  (let ((shape (dv-shape v)))
    (let loop ((q shape)
	       (indx '()))
      (cond ((null? q) 
	     (dv-set! v (reverse indx) 
			 (filler (reverse indx))))
	    (#t (let ((n (car q)))
		  (do ((j 0 (+ j 1))) ((= j n))
		    (loop (cdr q) (cons j indx)))))))))

;;; Sets values in v:
;;; (dv-set! v indx x): sets position `indx' in `v' with value `x'.  
;;; (dv-set! v filler): if filler is a procedure 
;;;             then sets position `indx' in `v' with value `(filler indx)'.  
;;; (dv-set! v x): sets all values of `v' to `x'
(define-public (dv-set! v . args)
  (cond ((= (length args) 2) (apply dv-set-w2 v args))
	(#t (let ((arg (car args)))
	      (cond ((procedure? arg) (dv-set-with-filler! v arg))
		    (#t (apply dv-set-w1 v args)))))))

;;; Sets w to the slice of `v' consisting of all elements
;;; such that index `indx' has a value in the set
;;; (start start+inc start+2*inc ...) that are lower than `end'. 
;;; range = indx start end [inc]
(define-public (dv-slice-range! v w range)
  (let ((shape (dv-shape w))
	(indx (car range))
	(start (cadr range))
	(end (caddr range))
	(inc 1))
    (if (>= indx (length shape)) (error "indx exceeds rank"))
    (if (>= (length range) 4) (set! inc (cadddr range)))
    (let ((v-shape shape)
	  (range-len (quotient (- end start) inc))
	  (filler (lambda (v-indx-vec) 
		    (let ((w-indx-vec v-indx-vec))
		      (list-set! w-indx-vec indx 
				 (+ start (* (list-ref v-indx-vec indx) inc)))
;		      (format #t "w-indx-vec ~A\n" w-indx-vec)
		      (dv-ref w w-indx-vec)))))
      (list-set! v-shape indx range-len)
      (apply dv-resize! v v-shape)
      (dv-set-with-filler! v filler))))

;;; range = indx start end [inc]
(define-public (dv-slice-indx! v w ivec indx)
  (let ((shape (dv-shape w)))
    (if (>= indx (length shape)) (error "indx exceeds rank"))
    (let ((v-shape shape)
	  (indx-len (dvint-size ivec)))
      (list-set! v-shape indx indx-len)
      (apply dv-resize! v v-shape)
      (dv-set-with-filler! v (lambda (v-indx-vec) 
				 (let ((w-indx-vec v-indx-vec))
				   (list-set! w-indx-vec indx 
					       (dvint-ref! ivec (list-ref v-indx-vec indx)))
				   (dv-ref w w-indx-vec)))))))

(define-public (dv-slice! v w . args)
  (cond ((pair? (car args)) (apply dv-slice-range! v w args))
	((dvint? (car args) (apply dv-slice-indx! v w args)))))


(define-public (dv-apply! v fun)
  (let ((n (dv-size v)))
    (do ((j 0 (+ j 1))) ((= j n)) 
      (let ((w (dv-ref v j)))
	(dv-set! v j (fun w))))))

(define-public (dv-add! v alpha)
  (dv-apply! v (lambda (y) (+ alpha y))))

(define-public (dv-rand! v)
  (dv-apply! v (lambda (y) (random:uniform))))

(define-public (dv-assoc v fun init)
  (let ((n (dv-size v)))
    (cond ((= n 0) init)
	  (#t (let ((result (dv-ref v 0)))
		(do ((j 1 (+ j 1))) ((= j n))
		  (set! result (fun (dv-ref v j) result)))
		result)))))

(define (dv-max-aux v comp)
  (dv-assoc v 
	       (lambda (x y) 
		 (cond ((comp x y) x) (#t y)))
	       0))

(define-public (dv-max v . args)
  (cond ((null? args) (dv-max-aux v >))
	(#t (dv-max-aux v (car args)))))
		
(define (dv-min-aux v comp)
  (dv-max v (lambda (x y) (not (comp x y)))))

(define-public (dv-min v . args)
  (cond ((null? args) (dv-max-aux v <))
	(#t (dv-min-aux v (car args)))))

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
  `(begin 
     (dv-method1 dvdbl ,fun)
     (dv-method1 dvint ,fun)))

(define-macro (dv-method-exp fun)
  `(begin 
     (dv-method1 dvdbl ,fun)
     (dv-method1 dvint ,fun)
     (export ,(dv-fun fun))))

(define-class <dvdbl> (<dvector>)
  (v #:accessor vec))

(define-method (initialize (v <dvdbl>) . args)
  (set! (vec v) (make-dvdbl)))

(dv-method resize-w!)
(dv-method set-w1)
(dv-method set-w2)

;(dv-method-exp clone!)

(define-method (dv-clone! (v <dvdbl>) (w <dvdbl>))
  (format #t "v ~A, w ~A\n" (vec v) (vec w))
  (dvdbl-clone! (vec v) (vec w)))

(define-method (dv-clone! (v <dvint>) (w <dvint>))
  (dvint-clone! (vec v) (vec w)))

(dv-method-exp push!)
(dv-method-exp size)
(dv-method-exp reshape!)
(dv-method-exp shape)
(dv-method-exp ref)
(dv-method-exp read!)
(dv-method-exp cat!)
(dv-method-exp dump)

(export <dvector> <dvdbl> vec 
	dv-resize! dv-set! dv-slice-range! 
	dv-slice-indx! dv-slice! 
	dv-apply! dv-add! dv-rand! dv-assoc dv-max-aux 
	dv-max dv-min dv-clone!)
