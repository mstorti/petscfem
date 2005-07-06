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
		(cond ((> x root) 
		       (list root (abb-load1 left x) right))
		      (else 
		       (list root left (abb-load1 right x))))))))

(define 

(format #t "abb ~A\n" (abb-load '() 3 4))
