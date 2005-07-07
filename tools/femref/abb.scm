(use-modules (ice-9 receive))

(define (abb-load abb1 . elems)
  (let loop ((inserted 0)
	     (abb abb1)
	     (q elems))
    (cond ((null? q) (values inserted abb))
	  (else 
	   (receive (ins1 abb-new) (abb-load1 abb (car q))
;		    (format #t "inserting ~A, inserted ~A\n"
;			    abb-new ins1)
		    (loop (+ ins1 inserted) abb-new (cdr q)))))))

(define (abb-load1 q x)
  (cond ((null? q) (values 1 (list x '() '())))
	(else (let ((root (car q))
		    (left (cadr q))
		    (right (caddr q)))
		(cond ((< x root) 
		       (receive (ins abb-new) (abb-load1 left x)
				(values ins (list root  abb-new right))))
		      ((> x root) 
		       (receive (ins abb-new) (abb-load1 right x)
				(values ins (list root left abb-new))))
		      (else (values 0 q)))))))

(define (abb-print abb1)
  (let loop ((abb abb1))
    (cond ((null? abb) (format #t " ."))
	  (else (let ((root (car abb))
		      (left (cadr abb))
		      (right (caddr abb)))
		  (cond ((and (null? left) (null? right)) (format #t " ~A" root))
			(else (format #t " (~A" root)
			      (loop left)
			      (loop right)
			      (format #t " )")))))))
  (newline))
  

(define (rand-list n m)
  (let loop ((q '())
	     (j 0))
    (cond ((= j n) q)
	  (else (loop (cons (random m) q) (+ j 1))))))

;; returns bool min max
(define (abb?-aux a)
  (cond ((null? a) (values #t #f #f))
	(else (let ((root (car a))
		    (left (cadr a))
		    (right (caddr a)))
		(receive (isl minl maxl) (abb?-aux left)
;			 (format #t "left ~A, isl ~A, minl ~A, maxl ~A\n"
;				 left isl minl maxl)
			 (cond ((not isl) (values #f #f #f))
			       ((and maxl (> maxl root)) 
;				(format #t "maxl ~A, root ~A\n" maxl root)
				(values #f #f #f))
			       (else 
				(receive (isr minr maxr) (abb?-aux right)
;					 (format #t "right ~A, isr ~A, minr ~A, maxr ~A\n"
;						 right isr minr maxr)
					 (cond ((not isr) (values #f #f #f))
					       ((and minr (< minr root)) (values #f #f #f))
					       (else 
;						(format #t "abb ~A, returns ~A ~A ~A\n" 
;							a #t minl maxr)
						(values #t 
							(cond (minl minl)
							      (else root))
							(cond (maxr maxr)
							      (else root)))))))))))))

(define (abb? a)
  (cond ((null? a) #t)
	(else (let ((root (car a))
		    (left (cadr a))
		    (right (caddr a)))
		(cond ((and (not (null? left)) (< (car left) root)) 
		       #f)
		      ((and (not (null? right)) (< (car right) root)) #f)
		      (else #t))))))

(define (tryme m n) 
  (let ((l (rand-list m n)))
    (format #t "inserting ~A\n" l)
    (receive (ins abb) (apply abb-load '() l)
	     (format #t "inserted ~A elements\n" ins)
	     (abb-print abb)))
  abb)

(define (abb-rand m n)
  (let ((l (rand-list m n)))
    (receive (ins abb) (apply abb-load '() l)
	     abb)))

;(tryme 5 4)
