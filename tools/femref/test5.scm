;;; $Id: test5.scm,v 1.12 2005/02/18 00:46:41 mstorti Exp $

(use-modules (srfi srfi-1))

;;---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
;;; Generate all partitions of 'n' We represente a partition as non
;;; increasing lists, for instance, the partitions of 3 are ((3) (2 1)
;;; (1 1 1)).
;;;
;;; A partial partition of `n' is a non increasing list whose sum is
;;; not greater than `n'.
;;;
;;; Given a partial partition `p' of `n' whose sum is `s' and whose
;;; smallest element if `m' all the partitions of `n' that start in `p'
;;; can be obtained by combaining `p' with all the partitions of `n-s'
;;; that have elements not greater than `m'.
;;;
;;; `(complete-1 part min-elem remain)' returns the list of partitions
;;; of `remain' in elements not greater than `min-elem'. Moreover,
;;; the partial partition `part' is appended to all these partitions. 
;;; So that the problem of finding the partitions of `n' starting with
;;; `part' can be solved as
;;; `(complete-1 part (apply min part) (- (n (apply sum part))))'. 
;;; And the global problem of finding all partitions of `n' can be solved as
;;; `(complete-1 '() n n)'. 
;;;
;;; Now, the partitions of `remain' in elements not greater than
;;; `min-elem' can be found by separating them in those who start in
;;; k=1, those who start with k=2, and so on, until `(min min-elem
;;; remain)'. The partitions of `remain' with elements not greater
;;; than `min-elem' that start with `k' can be obtained combining `k'
;;; with the partitions of `remain-k' in elements not greater than `k'.
;;; 
(define (partition n)
  ;; Completes one partial partition `part', with tails such that its
  ;; sum is `remain' and whose elements are not greater than
  ;; `min-elem'.
  (let complete-1 ((part '())
		   (min-elem n)
		   (remain n))
    (let ((kmax (min min-elem remain)))
      (cond ((zero? remain) (list part))
	    (#t 
	     ;; Loops over `k', cumulating lists in `completed-parts'. 
	     (let loop ((k 1)
			(completed-parts '()))
	       (cond ((> k kmax) completed-parts)
		     (#t 
		      (loop (+ k 1) 
			    (append
			     completed-parts 
			     ;; prepends `k' to part and pass to
			     ;; `complete-1' with rest `remain-k'. 
			     (complete-1 (cons k part) k (- remain k))))))))))))

(define (nparts2 n) (length (partition n)))

;;---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
;; Computes the number of partitions of `n'
(define (nparts n)
  ;; Auxiliary recursive function that computes the number of
  ;; partitions of `m' with elements not greater than `p'. This
  ;; algorithm is inefficient because it computes `(nparts-prtl m p)'
  ;; for the same pairt `(m p)' several times, instead of storing the
  ;; values. Some kind of lazy evaluation or simply storing the values
  ;; in a table (dynamic programming) would be much more efficient.
  (let nparts-prtl ((m n)
		    (p n))
    (cond ((= m 1) 1)
	  ((= m 0) 1)
	  ((= p 1) 1)
	  (#t 
	   (let ((k-max (min p m)))
	     (let loop ((q 0)
			(k 1))
	       (cond ((> k k-max) q)
		     (#t (loop (+ q (nparts-prtl (- m k) k)) (+ k 1))))))))))

;;---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
;; This is the more efficient algorithm, constructing the table.
;; In this form the algorithm is `n^3' (I think). 
(define (nparts3 n)
  (let* (
	 ;; This is the table that stores the previous values.
	 (table (make-array 0 n n))
	 ;; Auxiliary function that retrieves the value stored in the table
	 ;; for the pair `(n p)'
	 (np-ref (lambda (n p) 
		   (cond ((<= n 1) 1)
			 ((= p 1) 1)
			 (#t (array-ref table (- n 1) (- p 1))))))
	 ;; Auxiliary function that sets the entry in the table
	 ;; for the pair `(n p)'
	 (np-set! (lambda (val n p) 
		    (array-set! table val 
				(- n 1) (- p 1))))
	 ;; Computes one value in the table from previous values
	 (compute-1 (lambda (m p)
		      (let loop ((sub-parts 0)
				 (k 1))
			(cond ((> k (min p m)) sub-parts)
			      (#t (loop (+ sub-parts (np-ref (- m k) k)) (+ k 1))))))))
    ;; Loops over pairs `(m p)' increasing `p' first, until `(n n)' is
    ;; reached.
    (let loop ((m 1)
	       (p 1))
      ;; (format #t "m ~A, p ~A, table ~A\n" m p table)
      (cond ((> p m) (loop (+ m 1) 1))
	    ((> m n) (np-ref n n))
	    (#t 
	     (np-set! (compute-1 m p) m p)
	     (loop m (+ p 1)))))))

#!
(let loop ((n 1))
  (cond ((> n 15))
	(#t (format #t "(~A -> ~A ~A ~A)\n" 
		    n (nparts n) (nparts2 n) (nparts3 n))
	    (loop (+ n 1)))))

(let ((n 5))
      (format #t "(nparts3 ~A) ~A\n" n (nparts3 n)))

(let ((n 8))
      (format #t "(partition ~A) ~A\n" n (partition n)))
!#

