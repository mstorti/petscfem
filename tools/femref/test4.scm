;;; $Id: test4.scm,v 1.1 2005/02/14 04:27:47 mstorti Exp $

(define weight-rand
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
    (format #t "stat ~A\n" stat-l-f)))

(stat (weight-rand '((foo . 0.2) 
		     (bar . 0.1)
		     (baz . 0.3)
		     (foz . 0.4))) 50000)
