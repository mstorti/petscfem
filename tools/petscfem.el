;;; petscfem.el --- Mode for PETSc-FEM data files

;; Copyright (C) 2003  Free Software Foundation, Inc.

;; Author:  <mstorti@spider>
;; Keywords: languages, data

;; This file is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; This file is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with GNU Emacs; see the file COPYING.  If not, write to
;; the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
;; Boston, MA 02111-1307, USA.

;;; Commentary:

;; 

;;; Code:

(easy-mmode-defsyntax petscfem-mode-syntax-table
  '((?\# . "<")
   (?\n . ">"))
  "Syntax-table used in PETSc-FEM mode.")
;
(setq petscfem-font-lock-keywords 
      '(("^\\s-*\\(global_options\\|elemset\\|nodes\\|nodedata\\|fixa\\|fixa_amplitude\\|constraint\\|end_elemsets\\)\\>" . font-lock-type-face)
	("^\\w+\\>" . font-lock-variable-name-face)
	("^__\w*__\\>" . font-lock-keyword-face)
	("<:.*:>" . font-lock-warning-face)
	font-lock-keyword-face
	))
;
(define-derived-mode petscfem-mode
  text-mode "PETSc-FEM"
  "Major mode for editing PETSc-FEM data files
          \\{petscfem-mode-map}"
  (setq font-lock-defaults 
	`((petscfem-font-lock-keywords)
	  nil nil		  
	  ((?/ . "w") (?~ . "w") (?. . "w") (?- . "w") (?_ . "w")) nil))
;
;;   (setq font-lock-defaults 
;; 	`('(global_options elemset nodedata nodes)
;; 	  t nil		  
;; 	  ((?/ . "w") (?~ . "w") (?. . "w") (?- . "w") (?_ . "w")) nil))
;
;;   (setq font-lock-defaults
;; 	`((sh-font-lock-keywords
;; 	   sh-font-lock-keywords-1 sh-font-lock-keywords-2)
;; 	  nil nil
;; 	  ((?/ . "w") (?~ . "w") (?. . "w") (?- . "w") (?_ . "w")) nil
;; 	  (font-lock-syntactic-keywords
;; 	   ;; Copy so we can use destructive update in `sh-font-lock-heredoc'.
;; 	   . ,(copy-sequence sh-font-lock-syntactic-keywords))
;; 	  (font-lock-syntactic-face-function
;; 	   . sh-font-lock-syntactic-face-function)))  
  )
(setq petscfem-mode-hook 
      (lambda () 
	(setq info-lookup-mode 'petscfem-mode))
;
(load-library "info-look")
(info-lookup-maybe-add-help
 :mode 'petscfem-mode 
 :regexp "[_a-zA-Z0-9./+-]+"
 :doc-spec (list (list 
		  (concat "(" (getenv "PETSCFEM_DIR") 
			  "/doc/options.info)Options Index"))))
;
(defun my-Info-lookup-copy-keyword()
  (interactive)
  (save-excursion
    (goto-char (point-max))
    (search-backward "<<")
    (forward-char 2)
    (let ((beg (point)))
      (search-forward ">>")
      (forward-char -2)
      (princ (format "Yanked: \"%s\"" (buffer-substring beg (point))))
      (copy-region-as-kill beg (point)))))
;
(define-key Info-mode-map (kbd "<f4>") 'my-Info-lookup-copy-keyword)

(provide 'petscfem)
;;; petscfem.el ends here
