;;; $Id: utils.scm,v 1.1 2005/01/18 14:58:50 mstorti Exp $ 

(define-public (ddump s v)
  (format #t "\n\n~A: \n" s)
  (dvdbl-dump v))

(define (idump s v)
  (format #t "~A: \n" s)
  (dvint-dump v))

