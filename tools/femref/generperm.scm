(use-modules (ice-9 format))

(define (compose p1 p2)
  (let* ((n1 (vector-length p1))
	 (n2 (vector-length p2))
	 (q (make-vector n1))
	 (j 0))
    (if (not (= n1 n2)) (error "p1 and p2 should have equal length!\n"))
    (while (< j n1)
	   (let ((indx (vector-ref p1 j)))
	     (vector-set! q j (vector-ref p2 indx)))
	   (set! j (+ 1 j)))
    q))

(define (perm-id n)
  (let ((q (make-vector n))
	(j 0))
    (while (< j n)
	   (vector-set! q j j)
	   (set! j (+ 1 j)))
    q))

(define (prod-aux q G side)
  (let loop ((g G)
	     (qG '()))
    (cond ((null? g) qG)
	  (else (let ((pair-prod (cond ((eq? side 'left (compose q (car g))))
				       (else (compose (car g) q)))))
		  (if (member pair-prod G) 
		      (loop (cdr g) qG)
			(loop (cdr g) (cons pair-prod qG))))))))

(define (prod arg1 arg2)
  (cond ((list? arg2) (prod-aux arg1 arg2 'left))
	(else (prod-aux arg2 arg1 'right))))

;(define (generate q G)
;  (let loop ((GG 

;; (format #t "~A" (compose (vector 1 2 0) (vector 1 2 0)))

; (let ((n 10)
;       (j 0))
;   (while (< j n)
; 	 (format #t "~A\n" G)
; 	 (set! G (prod g G))
; 	 (set! j (+ 1 j))))

(define (lex-compare a b)
  (let ((na (vector-length a))
	(nb (vector-length b)))
    (if (not (= na nb)) (error "a and b should have same lenght\n"))
    (let loop ((j 0))
      (cond ((= j na) #f)
	    (else (let ((wa (vector-ref a j))
			(wb (vector-ref b j)))
		    (cond ((= wa wb) (loop (+ 1 j)))
			  (else (< wa wb)))))))))

(define (unique list less)
  (let ((sorted-list (sort list less)))
    (let loop ((q sorted-list)
	       (uniq-list '()))
      (cond ((null? q) (reverse uniq-list))
	    ((null? uniq-list) (loop (cdr q) (cons (car q) uniq-list)))
	    (else (cond ((or (less (car q) (car uniq-list))
			     (less (car uniq-list) (car q)))
			 (loop (cdr q) (cons (car q) uniq-list)))
			(else (loop (cdr q) uniq-list))))))))

(define (prod g G)
  (let ((gGg (append G 
		     (map (lambda (x) (compose g x)) G) 
		     (map (lambda (x) (compose x g)) G))))
    (unique gGg lex-compare)))

; (let ((a (vector 0 1 2 3 4))
;       (b (vector 0 1 1 3 4)))
;   (format #t "a: ~A, b: ~A, a<b: ~A\n" a b (lex-compare a b)))

; (let ((l '(1 2 3 1 2 3 4 5 6 4 5 6)))
;   (format #t "~A\n" (unique (reverse l) <)))

(define g (vector 1 2 3 4 5 0))
(define G (list g))

(define (generate-aux g G)
  (let loop ((p1 G)
	     (p2 (prod g G)))
	(cond ((= (length p1) (length p2)) p1)
	      (else (loop p2 (prod g p2))))))

(define (generate G)
  (if (zero? (length G)) (error "G can't be null\n"))
  (let loop ((gen-G (list (perm-id (vector-length (car G)))))
	     (g G))
    (cond ((null? g) gen-G)
	  (else (loop (generate-aux (car g) gen-G) (cdr g))))))

(define (check-gen gen-list s)
  (let* ((gen (generate gen-list))
	 (n (vector-length (car gen))))
    (format #t "total ~A perms: ~A, (of ~A! = ~A)\n" 
	    s (length gen) n (factorial n))
    (map (lambda (x) (format #t "~A\n" x)) gen)))

(define (factorial n)
  (let loop ((fact 1)
	     (count 1))
    (cond ((> count n) fact)
	  (else (loop (* fact count) (+ 1 count))))))

; (format #t "(factorial 5): ~A\n" (factorial 5))

(check-gen (list (vector 1 2 0)) "tri")

(check-gen (list (vector 1 2 0)
		 (vector 1 0 2)) "unoriented tri")

(check-gen (list (vector 1 2 3 0)) "quad")

(check-gen (list (vector 1 2 3 0)
		 (vector 0 3 2 1)) "unoriented quad")

(check-gen (list (vector 1 3 2 0) 
		 (vector 1 2 0 3)) "tetra")

(check-gen (list (vector 1 3 2 0) 
		 (vector 1 2 0 3)
		 (vector 1 0 2 3)) "unoriented tetra")

(check-gen (list (vector 1 2 3 0 5 6 7 0)
		 (vector 1 5 6 2 0 4 7 3)) "hexa")

(check-gen (list (vector 1 2 3 0 5 6 7 0)
		 (vector 1 5 6 2 0 4 7 3)
		 (vector 1 0 3 2 5 4 7 6)) "unoriented hexa")

(check-gen (list (vector 4 3 5 1 0 2)
		 (vector 1 2 0 4 5 3)) "oriented prism")

(check-gen (list (vector 4 3 5 1 0 2) ; quad face rot
		 (vector 1 2 0 4 5 3) ; tri face rot
		 (vector 1 0 2 4 3 5)) "unoriented prism")
