;;; $Id: dvector.scm,v 1.24 2005/10/17 02:51:37 mstorti Exp $
(define-module (dvector))
(use-modules (oop goops))
(use-modules (ice-9 optargs))

(load-extension 
 (search-path %load-path "libfemref.so") 
 "dvint_init")
(load-extension 
 (search-path %load-path "libfemref.so") 
 "dvdbl_init")

(define dv-version "$Id: dvector.scm,v 1.24 2005/10/17 02:51:37 mstorti Exp $")

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
;;; range = start end [inc]
;;; If `start' is `#f' then 0 is used.
;;; If `end' is `#f' then the size of `indx' index is used. 
;;; Example: (dv-slice-range! v w '(0 5 2) 1)
;;; extracts all elements such that index 1 takes the values
;;; 0 2 4 
(define-public (dv-slice-range! v w range indx)
  (let ((shape (dv-shape w))
	(start (car range))
	(end (cadr range))
	(inc 1))
    (if (not start) (set! start 0))
    (if (not end) (set! end (list-ref shape indx)))
    (if (>= indx (length shape)) (error "indx exceeds rank"))
    (if (>= (length range) 4) (set! inc (cadddr range)))
    (let ((v-shape shape)
	  (range-len (+ 1 (quotient (- end start 1) inc)))
	  (filler (lambda (v-indx-vec) 
		    (let ((w-indx-vec v-indx-vec))
		      (list-set! w-indx-vec indx 
				 (+ start (* (list-ref v-indx-vec indx) inc)))
		      (dv-ref w w-indx-vec)))))
      (list-set! v-shape indx range-len)
      (apply dv-resize! v v-shape)
      (dv-set-with-filler! v filler))))

;;; Sets `w' to the slice of `v' consisting of all
;;; elements for which index takes the values given in `ivec'.
;;; ivec is a `<dvint>' vector. 
(define-public (dv-slice-indx! v w ivec indx)
  (if (not (is-a? ivec <dvint>)) 
      (error "ivec is not an <dvint> vector. ivec= ~A" ivec))
  (let ((shape (dv-shape w)))
    (if (>= indx (length shape)) (error "indx exceeds rank"))
    (let ((v-shape shape)
	  (indx-len (dv-size ivec))
	  (filler (lambda (v-indx-vec) 
		    (let ((w-indx-vec v-indx-vec))
		      (list-set! w-indx-vec indx 
				 (dv-ref ivec (list-ref v-indx-vec indx)))
		      (dv-ref w w-indx-vec)))))
      (list-set! v-shape indx indx-len)
      (apply dv-resize! v v-shape)
      (dv-set-with-filler! v filler))))

;;; Auxiliary function that calls either the
;;; `-range!' or `index!' versions. 
(define (dv-slice1! v w range indx)
  (cond ((is-a? range <dvint>)
	 (dv-slice-indx! v w range indx))
	(else
	 (dv-slice-range! v w range indx))))

;;; Wrapper that iterates over pair of arguments appliying
;;; successively `dv-slice1!'.
;;; Example:
;;;    (dv-slice! v w '(0 5 2) 0 '(0 3) 1)
(define-public (dv-slice! v w . args)
  (cond ((null? args) (dv-clone! v w))
	((= (length args) 1) (error "[1] not enough args"))
	((= (length args) 2) (apply dv-slice1! v w args))
	(else 
	 ;;; Each pair of values in `args' defines a slicing.
	 ;;; For instance, if we write
	 ;;;    (dv-slice! v w r1 i1 r2 i2 r3 i3 r4 i4)
	 ;;; then this should be expanded to
	 ;;;    (dv-slice! vaux1 w r1 i1)
	 ;;;    (dv-slice! vaux2 vaux1 r2 i2)
	 ;;;    (dv-slice! vaux3 vaux2 r3 i3)
	 ;;;    (dv-slice! v vaux3 r4 i4)
	 ;;; But in fact we need only two auxiliary vectors
	 ;;; `vaux1', `vaux2' alternating in the following way
	 ;;;    (dv-slice! vaux1 w r1 i1)
	 ;;;    (dv-slice! vaux2 vaux1 r2 i2)
	 ;;;    (dv-slice! vaux1 vaux2 r3 i3)
	 ;;;    (dv-slice! v vaux1 r4 i4)
	 ;;; The variables `v1' and `v2' below takes the alternate
	 ;;; values (vaux1,vaux2), (vaux2,vaux1) and so on.
	 ;;; When the numbr of args reaches 2, then we directly
	 ;;; slice onto `v'. 
	 (let ((vaux1 (make (class-of w)))
	       (vaux2 (make (class-of w))))
	   (dv-slice1! vaux1 w (car args) (cadr args))
	   (let loop ((q (cddr args))
		      (v1 vaux1)
		      (v2 vaux2))
	     (cond ((= (length q) 1) (error "[2] not enough args"))
		   ((= (length q) 2) 
		    (dv-slice1! v v1 (car q) (cadr q)))
		   (else 
;		    (format #t "q ~A\n" q)
		    (dv-slice1! v2 v1 (car q) (cadr q))
		    (loop v1 v2 (cddr q)))))))))

;;; Applys a function to each element in `v'.
;;; Example: (dv-apply! v (lambda(x) (* 2 x)))
(define-public (dv-apply! v fun)
  (let ((n (dv-size v)))
    (do ((j 0 (+ j 1))) ((= j n)) 
      (let ((w (dv-ref v j)))
	(dv-set! v j (fun w))))))

;;; Adds a fixed number to each element in `v'. 
;;; Example: (dv-add! v 1.3)
(define-public (dv-add! v alpha)
  (dv-apply! v (lambda (y) (+ alpha y))))

;;; usage: (dv-reduce v fun)
;;;        (dv-reduce v fun knil)
;;; Reduces all elements with function `fun'.
;;; `fun' should be a function that is invariant under
;;; permutations of its arguments, and furthermore
;;; we assume that there is a value `knil' such that
;;; r = (f v0 v1 ...) is equivalent to
;;; r = knil; r=(fun v_0 r), r=(fun v_1 r), ...
;;; For instance, for `fun=+' we have `knil=0', and also
;;; (fun=max,knil=-infty), (fun=min,knil=+infty),
;;; (fun=max_abs,knil=0), (fun=min_abs,knil=+infty),
;;; (fun=sum_square,knil=0), (fun=sum_abs,knil=0), 
;;; (fun=*,knil=1). 
;;; In the first form, i.e. if `knil' is not provided,
;;; then `fun' should accept 0, 1 or two arguments. If
;;; the vector has at less two arguments, then `fun' is
;;; guaranteed to be called only with two arguments.
;;; If it is called with 1 or 0 arguments, then `fun' should
;;; accept the same arity.
;;; In the second form, i.e. if `knil' is provided,
;;; then  `fun' is
;;; guaranteed to be called only with two arguments.
;;; `knil' is the neutral value for `fun', i.e. (fun args)
;;; should be invariant if we insert `knil' in each position of
;;; `args'. 
;;; Example: (dv-reduce v (lambda(x y) (+ x y)) 0)
;;; or:
;;;     (dv-reduce v (lambda args 
;;; 	       (cond ((zero? (length args)) 0)
;;; 		     ((= (length args) 1) (car args))
;;; 		     (else (apply + args)))))
(define-public (dv-reduce v fun . args)
  (cond ((not (null? args))
	 (apply dv-reduce-nil v fun args))
	(else
	 (dv-reduce-non-nil v fun))))

;;; This is for the case of non having a `nil' object.
;;; The function then can be called with zero one or
;;; two arguments.
(define (dv-reduce-non-nil v fun)
  (let ((n (dv-size v)))
    (cond ((= n 0) (fun))
	  ((= n 1) (fun (dv-ref v 0)))
	  (else
	   (let ((n (dv-size v)))
	     (let loop ((j 2)
			(result (fun (dv-ref v 0) (dv-ref v 1))))
	       (cond ((= j n) result)
		     (else (loop (+ j 1) (fun result (dv-ref v j)))))))))))

;;; This is when having a nil object. 	    
(define (dv-reduce-nil v fun knil)
  (let ((n (dv-size v)))
    (let loop ((j 0)
	       (result knil))
      (cond ((= j n) result)
	    (else (loop (+ j 1) (fun result (dv-ref v j))))))))

;;; Prints vector to a file or standard output. 	    
;;; usage:
;;;    (dv-dump v . file message items)
(define-public (dv-dump v . args)
  (let-keywords args #f ((file  #f)
 			 (rowsz #f))
;;;   (let-optional args ((file  #f)
;;; 		      (rowsz #f))
     (format #t "file ~A, rowsz ~A, args ~A\n" file rowsz args)
;      (if (null? args) 
; 	(begin 
; 	  (dv-dump1 v file rowsz))
; 	(begin 
; 	  (if (= (length args) 1)
; 	      (format #t "~A: \n" (car args))
; 	      (apply format #t args))
; 	  (dv-dump1 v file rowsz)
; 	  (newline)))
))

(define (dv-max-aux v comp)
  (dv-reduce v (lambda (x y) 
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

(define-class <dvint> (<dvector>)
  (v #:accessor vec))

(define-method (initialize (v <dvint>) . args)
  (set! (vec v) (make-dvint)))

(define-method (dv-clone! (v <dvdbl>) (w <dvdbl>))
  (dvdbl-clone! (vec v) (vec w)))

(define-method (dv-clone! (v <dvint>) (w <dvint>))
  (dvint-clone! (vec v) (vec w)))

;;; Fills with uniform random numbers (uniform
;;; distribution in [0,1])
(define-method (dv-rand! (v <dvdbl>))
  (dv-apply! v (lambda (x) (random:uniform))))

;;; Fills with uniform random numbers in the
;;; range [0,n)
(define-method (dv-rand! (v <dvint>) n)
  (dv-apply! v (lambda (x) (random n))))

(dv-method resize-w!)
(dv-method set-w1)
(dv-method set-w2)
(dv-method dump1)

(dv-method-exp push!)
(dv-method-exp size)
(dv-method-exp reshape!)
(dv-method-exp shape)
(dv-method-exp ref)
(dv-method-exp read!)
(dv-method-exp cat!)

(export <dvector> <dvdbl> <dvint> vec 
	dv-resize! dv-set! dv-slice-range! 
	dv-slice-indx! dv-slice! 
	dv-apply! dv-add! dv-rand! dv-reduce dv-max-aux 
	dv-max dv-min dv-clone! dv-slice!
	dv-version dv-rand! make-dvint make-dvdbl)
