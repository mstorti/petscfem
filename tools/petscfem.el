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

;; Basic support for editing PETSc-FEM data files with Emacs. Supports
;; some basic colorization (complex ePerl constructs nay fool it), and
;; also some code to automatically insert PETSc-FEM keywords in the
;; data file. 

;;; Code:

(easy-mmode-defsyntax petscfem-mode-syntax-table
  '((?\# . "<")
   (?\n . ">"))
  "Syntax-table used in PETSc-FEM mode.")
;
(defvar petscfem-font-lock-keywords 
  '(("^\\s-*\\(global_options\\|elemset\\|nodes\\|nodedata\\|fixa\\|fixa_amplitude\\|constraint\\|end_elemsets\\)\\>" . font-lock-type-face)
    ("^\\w+\\>" . font-lock-variable-name-face)
    ("^__\w*__\\>" . font-lock-keyword-face)
    ("<:.*?:>" . font-lock-warning-face)
    ("^>>\\(if\\|elsif\\|endif\\|else\\)\\>" . font-lock-warning-face)
;;    ("\#if" . font-lock-warning-face)
    font-lock-keyword-face
    )
  "Basic colorization scheme for PETSc-FEM. ")
;;
(defvar petscfem-dir
  (getenv "PETSCFEM_DIR")
  "The root directory of the PETSc-FEM installation.")
;;
(defvar petscfem-info-file
  (when (boundp 'petscfem-dir)
    (concat petscfem-dir "/doc/options.info"))
  "The PETSc-FEM options info file for using with info-lookup.")
;;
(define-derived-mode petscfem-mode
  text-mode "PETSc-FEM"
  "Major mode for editing PETSc-FEM data files
          \\{petscfem-mode-map}."
  (setq font-lock-defaults 
	`((petscfem-font-lock-keywords)
	  nil nil		  
	  ((?/ . "w") (?~ . "w") (?. . "w") (?- . "w") (?_ . "w")) nil))
 )
;;
(defvar petscfem-mode-hook nil 
  "Normal hook run when entering PETSc-FEM.")
;;
(defun my-Info-lookup-copy-keyword()
  (interactive)
  (save-excursion
    (when (not (search-forward ">>"))
      (error "Probably not in a PETSc-FEM options info buffer!!"))
      (forward-char -2)
    (let ((end (point)))
      (when (not (search-backward "<<"))
	(error "Probably not in a PETSc-FEM options info buffer!!"))
      (forward-char 2)
      (princ (format "Yanked: \"%s\"" (buffer-substring (point) end)))
      (copy-region-as-kill (point) end))))
;
(defun my-Info-bury-and-kill()
  (interactive)
  (Info-exit)
  (delete-window))
;
(provide 'petscfem)
;;; petscfem.el ends here
