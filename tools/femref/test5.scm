;;; $Id: test5.scm,v 1.4 2005/02/16 01:38:30 mstorti Exp $

(use-modules (srfi srfi-1))

(define (sum l) (fold (lambda (a b) (+ a b)) 0 l))

;;; ================================================================
;;; Generate all partitions of 'n'
(define (partition n)
  ;;; Completes one partition
  (define (complete-1 part min-elem remain)
    (format #t "part ~A, min-elem ~A, remain ~A\n"
	    part min-elem remain)
    (let ((kmax (min min-elem remain)))
      (cond ((zero? remain) (list part))
	    (#t 
	     (let loop ((k 1)
			(completed-parts '()))
	       (format #t "k ~A, completed-parts ~A\n"
		       k completed-parts)
	       (cond ((> k kmax) 
		      (format #t 
			      "loop returning completed-parts ~A\n" 
			      completed-parts) 
		      completed-parts)
		     (#t 
		      (loop (+ k 1) 
			    (append
			     completed-parts 
			     (complete-1 (cons k part) k (- remain k)))))))))))
  (complete-1 '() n n))


(let ((n 6))
  (format #t "part ~A: ~A" n (partition n)))
