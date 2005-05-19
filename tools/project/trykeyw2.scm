;;; $Id: trykeyw2.scm,v 1.1 2005/05/19 02:26:01 mstorti Exp $
(use-modules (ice-9 format))
(use-modules (ice-9 receive))

(define (parse-args kws args)
  (let ((kw-vals (list-copy kws)))
    (let loop ((q args)
	       (clean-args '()))
      (cond ((null? q) (values kw-vals clean-args))
	    (else 
	     (let ((key (car q)))
	       (cond ((not (symbol? key)) (loop (cdr q) (cons key clean-args)))
		     (else (let ((qk (assq key kw-vals)))
			     (cond (qk (if (null? (cdr q)) 
					   (error (format #f "keyw with no arg ~A" key)))
				       (set-cdr! qk (cadr q))
				       (loop (cddr q) clean-args))
				   (else (loop (cdr q) (cons key clean-args)))))))))))))

(receive 
 (kw-vals args)
 (parse-args '((a . 1) (b . 5)) '(34 45 3 b 7 w a))
 (format #t "kws ~A, args ~A\n" kw-vals args))
