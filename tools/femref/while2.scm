;; The inner `do' loop avoids re-establishing a catch every iteration,
;; that's only necessary if continue is actually used.  A new key is
;; generated every time, so break and continue apply to their originating
;; `while' even when recursing.  `while-helper' is an easy way to keep the
;; `key' binding away from the cond and body code.
;;
(define-macro (while cond . body)
  (define (while-helper proc)
    (do ((key 'my-secret-while-key))
	((catch key
		(lambda ()
		  (proc (lambda () (throw key #t))
			(lambda () (throw key #f))))
		(lambda (key arg) arg)))))
  `(,while-helper (,lambda (break continue)
		    (,do ()
			((,not ,cond))
		      ,@body)
		    #t)))
