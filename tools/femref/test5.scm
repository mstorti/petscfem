;;; $Id: test5.scm,v 1.8 2005/02/17 00:39:09 mstorti Exp $

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
  (define (nparts-aux n p)
;;;  (format #t "n ~A, p ~A\n" n p)
    (cond ((= n 1) 1)
	((= n 0) 1)
	((= p 1) 1)
	(#t 
	 (let ((k-max (min p n)))
	   (let loop ((m 0)
		      (k 1))
;;;	     (format #t "m ~A, k ~A\n" m k)
	     (cond ((> k k-max) m)
		   (#t (loop (+ m (nparts-aux (- n k) k)) (+ k 1)))))))))
  (nparts-aux n n))

(define (nparts2 n) (length (partition n)))

#!
(let* ((n 7)
       (parts (partition n)))
  (format #t "part ~A: total ~A partitions.\n" n (length parts))
;  (format #t "partitions: ~A\n" parts)
  (format #t "(nparts n) -> ~A\n" (nparts n n)))
!#

#!
(let loop ((n 1)
	   (p 1))
  (cond ((> p n) (format #t "\n") (loop (+ n 1) 1))
	((> n 5))
	(#t (format #t "(~A ~A -> ~A ~A) " n p (nparts n p))  
	    (loop n (+ p 1)))))
!#

(define (nparts3 n)
  (let ((table (make-array 0 n n)))
    (let loop ((m 1)
	       (p 1))
      (cond ((> p n) (loop (+ n 1) 1))
	    ((> m n) (array-ref table n n))
	    (#t 
	     (array-set! table
	      (let loop2 ((sub-parts 0)
			  (k 1))
		(cond ((> k (min p n)) sub-parts)
		      (#t (loop2 (+ sub-parts (array-ref table (- m k) k))))))
			  m p))))))

#!
(let loop ((n 1))
  (cond ((> n 15))
	(#t (format #t "(~A -> ~A ~A)\n" n (nparts n) (nparts2 n))
	    (loop (+ n 1)))))
!#

(let ((n 8))
      (format #t "(nparts3 ~A) ~A\n" n (nparts2 n)))
