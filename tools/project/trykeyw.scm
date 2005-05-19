;;; $Id: trykeyw.scm,v 1.3 2005/05/19 11:56:49 mstorti Exp $
(use-modules (ice-9 format))
(use-modules (ice-9 optargs))

;;; ttn/optargs-kw-utils.scm --- make using (ice-9 optargs-kw) easier
;; Copyright (C) 2003-2005 Thien-Thi Nguyen
(define (remove-keys ls)
  (let loop ((ls ls) (acc '()))
    (if (null? ls)
        (reverse! acc)
        (let ((kw? (keyword? (car ls))))
          (loop ((if kw? cddr cdr) ls)
                (if kw? acc (cons (car ls) acc)))))))

(define* (f x y #:key a b . rest)
  (set! rest (remove-keys rest))
  (format #t "x ~A, y ~A, a ~A, b ~A, rest ~A\n"
	  x y a b rest))

(f 1 2 3)
(f 1 2 3 #:a 4)
(f 1 2 #:a 4 3)

(define* (g #:key a b)
	   (format #t "a ~A, b ~A\n" a b))

(g #:a 1 #:b 2)
(g #:b 1 #:a 2)


