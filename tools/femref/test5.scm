;;; $Id: test5.scm,v 1.2 2005/02/16 01:22:56 mstorti Exp $

(use-modules (srfi srfi-1))

(define (sum l) (fold (lambda (a b) (+ a b)) 0 l))

;;; ================================================================
;;; Generate all partitions of 'n'
(define (partition n)
  ;;; Completes one partition
  (define (complete-1 part min-elem remain)
    (let ((kmax (min min-elem remain)))
      (cond ((zero? remain) (list part))
	    (#t 
	     (let loop ((k 1)
			(completed-parts '()))
	       (cond ((= k kmax) completed-parts)
		     (#t 
		      (loop (+ k 1) 
			    (map (lambda (q) 
				   (append q part))
				 (append
				  completed-parts 
				  (complete-1 (cons k part) k (- n k))))))))))))
  (complete-1 '() n n))

#!
;;; ================================================================
;;; Generate all partitions of 'n'
(define (partition2 n)
  (let (part-list '())		; List of partial partitions
    (let loop ((q part-list)
	       (completed-lists '()))
      (cond ((null? q) completed-lists)
	    (#t (append 
		 completed-lists 
		 ;; Given a `partial list' returns the partitions
		 ;; originated in from it. 
		 (let loop2 ((partial-list (car q))
			     (result '()))
		   (let ((min 


  (let loop ((part-list '(()))
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
!#

(let ((n 4))
  (format #t "part ~A: ~A" n (partition n)))
