;;; $Id: test2.scm,v 1.1 2005/01/18 20:38:17 mstorti Exp $

(define-macro (def10 s) `(define ,s 10))

(define q 20)
(def10 q)
(format #t "q ~A\n" q)

(define qq-rr 0)
;(format #t "~A\n" (string->symbol (string-append (symbol->string 'qq) "-" "rr")))

; (define-macro (my-def a b c) `(define ,(string->symbol (string-append "qq" "rr")) c))
; (my-def qq rr 20)
; (format "qq-rr ~\n" qq-rr)

(define-macro (my-def a)
  `(define ,(string->symbol (string-append (symbol->string a) "-rr")) 40))
(my-def qq)
(format #t "qq-rr ~A\n" qq-rr)
