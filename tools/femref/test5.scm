;;; $Id: test5.scm,v 1.6 2005/02/16 23:41:07 mstorti Exp $

(use-modules (srfi srfi-1))

;;; ================================================================
;;; Generate all partitions of 'n'
(define (partition n)
  ;;; Completes one partition
  (let complete-1 ((part '())
		   (min-elem n)
		   (remain n))
;     (format #t "part ~A, min-elem ~A, remain ~A\n"
; 	    part min-elem remain)
    (let ((kmax (min min-elem remain)))
      (cond ((zero? remain) (list part))
	    (#t 
	     (let loop ((k 1)
			(completed-parts '()))
; 	       (format #t "k ~A, completed-parts ~A\n"
; 		       k completed-parts)
	       (cond ((> k kmax) 
; 		      (format #t 
; 			      "loop returning completed-parts ~A\n" 
; 			      completed-parts) 
		      completed-parts)
		     (#t 
		      (loop (+ k 1) 
			    (append
			     completed-parts 
			     (complete-1 (cons k part) k (- remain k))))))))))))

(define (nparts n)
  (cond ((= n 1) 1)
	(#t (let loop ((m 0)
		       (k 1))
	      (cond ((> k n) m)
		    (#t (loop (+ m (* k (nparts (- n k)))) (+ k 1))))))))


	

(let* ((n 4)
       (parts (partition n)))
  (format #t "part ~A: total ~A partitions.\n" n (length parts))
  (format #t "partitions: ~A\n" parts)
  (format #t "(nparts n) -> ~A\n" (nparts n)))
