(define (abb-load abb1 . elems)
  (let loop ((abb abb1)
	     (q elems))
    (cond ((null? q) abb)
	  (else (loop (abb-load1 abb (car q)) (cdr q))))))

(define (abb-load1 q x)
  (cond ((null? q) (list x '() '()))
	(else (let ((root (car q))
		    (left (cadr q))
		    (right (caddr q)))
		(cond ((< x root) 
		       (list root (abb-load1 left x) right))
		      (else 
		       (list root left (abb-load1 right x))))))))

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

;(define abb (abb-load '() 3 4 5 6 7 8))

(define (tryme m n) 
  (abb-print (apply abb-load '() (rand-list m n))))


