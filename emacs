; -*- mode: lisp -*-
(setq next-line-add-newlines nil) ; para que no apendice lineas
(defun fastmat2-print()
  (interactive)
  (let ((ident (buffer-substring 
		(+ (re-search-backward "[^0-9a-zA-Z_]") 1)
		(+ (re-search-forward 
		    "[^0-9a-zA-Z_]" (point-max-marker) t 2) -1)
		)))
    (switch-to-buffer "*gud*")
    (insert "p " ident ".print(\"\"),1")
    (comint-send-input)
    )
  )


(add-hook 'c++-mode-hook '(lambda () 
      (abbrev-mode 1)
      (auto-fill-mode)
;      (local-set-key "\C-cc" 'compile)
      (local-set-key "\C-ck" 'kill-compilation)
      (local-set-key [f10] 'fastmat2-print)
      (setq current-pos-regexp "^ *// Current line ===========")
       ))

;;
;; Reemplazar identificadores
;;
(defun query-replace-identifier()
  (interactive)
  (query-replace-regexp (concat "\\([^a-zA-Z]\\)"
				(read-from-minibuffer "replace identifier: ") 
				"\\([^a-zA-Z0-9]\\)")
			(concat "\\1" (read-from-minibuffer "with: ") 
				 "\\2")
			)
)
;
; This way it works for any kind of machine name
;
(defalias 'notas-insert-date (read-kbd-macro
			      "%=== <<date>> C-e SPC C-u ESC !uname SPC -n RET <down> C-k \
<up> C-e C-r SPC RET <right> C-u ESC !whoami RET <down> C-a DEL @ C-e <down>"))

;;
;;
(defalias 'pqpq
  (read-kbd-macro "C-x n n <home> <home> C-SPC <end> C-x C-x C-M-% , SPC * C-q LFD RET C-q LFD RET ! <home> ESC % , RET C-q LFD RET ! <end> C-x n w C-e ;"))

;;  (read-kbd-macro "C-x n n <home> C-M-% , SPC * C-q LFD RET ! RET <home> ESC % , RET C-q LFD RET ! RET <home> C-x n w"))
