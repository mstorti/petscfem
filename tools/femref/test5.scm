;;; $Id: test5.scm,v 1.1 2005/02/16 00:25:58 mstorti Exp $

;;; Generate all partitions of 'n'
(define (partition n)
  (let loop ((partial-list '())
	     (min-elem n)
	     (remain n))
    (cond ((= remain 0) partial-list)
	  (#t 
	   (format #t "pl ~A, min-elem ~A, remain ~A\n" 
		   partial-list min-elem remain)
	   (let loop2 ((completed-list '())
		       (k 1))
	     (format #t "cl ~A, k ~A\n" completed-list k)
	     (cond ((> k (max remain min-elem)) completed-list)
		   (#t (loop2 
			(append completed-list 
				(map (lambda (q) 
				       (format #t "q ~A, pl ~A\n"
					      q partial-list)
				       (append q partial-list))
				     (loop (cons k partial-list) k (- remain k))))
			(+ k 1)))))))))

(partition 4)
