;;; $Id: test4.scm,v 1.2 2005/02/14 22:16:30 mstorti Exp $

(define weight-rand-l
  (lambda (w-alist-a)
    (let ((w-alist (list-copy w-alist-a)))
      (lambda ()
	(let ((x (random:uniform)))
	  (let loop ((q w-alist)
		     (cum-w 0))
	    (let ((cum-w-new (+ cum-w (cdar q))))
	      (cond ((null? q) (error "not normalized weights"))
		    ((< x cum-w-new) (caar q))
		    (#t (loop (cdr q) cum-w-new))))))))))

(define (weight-rand-f f x1 x2 n)
  (let ((h (/ (- x2 x1) n)))
    (let loop ((w-alist '())
	       (x x1)
	       (cum-w 0))
      (let ((fx (f x)))
	(cond ((> x x2)
	       (weight-rand-l 
		(reverse 
		 (map (lambda (x) (cons (car x) (/ (cdr x) cum-w)))
		      w-alist))))
	      (#t (loop (cons (cons x fx) w-alist) (+ x h) (+ cum-w fx))))))))
	     
(define (weight-rand . args)
  (cond ((pair? (car args)) (weight-rand-l (car args)))
	((procedure? (car args))
	 (apply weight-rand-f args))
	(#t (error "invalid arguments\n"))))

(define (stat rand-gen n)
  (let ((stat-l-f
	 (let loop ((k n)
		     (stat-l '()))
	   (let* ((x (rand-gen))
		  (q (assoc x stat-l)))
	     (cond ((zero? k) stat-l)
		   ((not q) (loop (- k 1) (cons (cons x 1) stat-l)))
		   (#t (set-cdr! q (+ 1 (cdr q)))
		       (loop (- k 1) stat-l)))))))
    (map (lambda (q) (set-cdr! q (/ (cdr q) n))) stat-l-f)
    stat-l-f))

#!
(format #t "stat: ~A\n"
	(stat (weight-rand '((foo . 0.2) 
			     (bar . 0.1)
			     (baz . 0.3)
			     (foz . 0.4))) 5))

(map (lambda (x) (format #t "~A ~A\n" (car x) (cdr x)))
     (sort (stat (weight-rand (lambda (x) x) 0 1 50) 2000)
	   (lambda (a b) (< (car a) (car b)))))
!#

(format #t "Stat for gaussian prob fun\n")
(map (lambda (x) (format #t "~A ~A\n" (car x) (cdr x)))
     (sort (stat (weight-rand (lambda (x) (exp (- (* x x)))) -3 3 50) 20000)
	   (lambda (a b) (< (car a) (car b)))))

