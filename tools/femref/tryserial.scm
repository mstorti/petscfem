;;; $Id: tryserial.scm,v 1.1 2005/02/07 13:13:27 mstorti Exp $ 
;;; Example of serialization of objects

(define u '(5 4 3 2 #u(1 2 3 4 5)))
(define up (format #f "~A" u))
(format #t "printed representation: ~A\n" up)
(define uu (with-input-from-string up read))
(format #t "(eq? u uu) ~A\n" (eq? u uu))
(format #t "(eqv? u uu) ~A\n" (eqv? u uu))
(format #t "(equal? u uu) ~A\n" (equal? u uu))

