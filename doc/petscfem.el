;;; petscfem-mode.el --- text mode, and its idiosyncratic commands

;; Copyright (C) 1985, 1992, 1994 Free Software Foundation, Inc.

;; Maintainer: FSF
;; Keywords: wp

;; This file is part of GNU Emacs.

;; GNU Emacs is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; GNU Emacs is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with GNU Emacs; see the file COPYING.  If not, write to the
;; Free Software Foundation, Inc., 59 Temple Place - Suite 330,
;; Boston, MA 02111-1307, USA.

;;; Commentary:

;;; Code:

(defcustom petscfem-mode-hook nil
  "Hook run when entering PETSc-FEM."
  :type 'hook
  :group 'data)

(defvar petscfem-mode-variant nil
  "Non-nil if this buffer's major mode is a variant of PETSc-FEM mode.")

(defvar petscfem-mode-syntax-table nil
  "Syntax table used while in PETSc-FEM mode.")

(defvar petscfem-mode-abbrev-table nil
  "Abbrev table used while in petscfem mode.")

(define-abbrev-table 'petscfem-mode-abbrev-table ())

(if petscfem-mode-syntax-table
    ()
  (setq petscfem-mode-syntax-table (copy-syntax-table sh-mode-syntax-table)))

(defvar petscfem-mode-map nil
  "Keymap for PETSc-FEM mode.")

(if petscfem-mode-map
    ()
  (setq petscfem-mode-map (make-sparse-keymap)))


(defun petscfem-mode ()
  "Major mode for editing PETSc-FEM data files. 
\\{petscfem-mode-map}
Turning on PETSc-FEM mode runs the normal hook `petscfem-mode-hook'."
  (interactive)
  (kill-all-local-variables)
  (use-local-map petscfem-mode-map)
  (setq local-abbrev-table petscfem-mode-abbrev-table)
  (set-syntax-table petscfem-mode-syntax-table)
  (make-local-variable 'paragraph-start)
  (setq paragraph-start (concat page-delimiter "\\|[ \t]*$"))
  (if (eq ?^ (aref paragraph-start 0))
      (setq paragraph-start (substring paragraph-start 1)))
  (make-local-variable 'paragraph-separate)
  (setq paragraph-separate paragraph-start)
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'indent-relative-maybe)
  (setq mode-name "PETSc-FEM")
  (setq major-mode 'petscfem-mode)
  (run-hooks 'petscfem-mode-hook))

(defun paragraph-indent-petscfem-mode ()
  "Major mode for editing text, with leading spaces starting a paragraph.
In this mode, you do not need blank lines between paragraphs
when the first line of the following paragraph starts with whitespace.
`paragraph-indent-minor-mode' provides a similar facility as a minor mode.
Special commands:
\\{petscfem-mode-map}
Turning on Paragraph-Indent Text mode runs the normal hooks
`petscfem-mode-hook' and `paragraph-indent-petscfem-mode-hook'."
  (interactive)
  (kill-all-local-variables)
  (use-local-map petscfem-mode-map)
  (setq mode-name "Parindent")
  (setq major-mode 'paragraph-indent-petscfem-mode)
  (setq local-abbrev-table petscfem-mode-abbrev-table)
  (set-syntax-table petscfem-mode-syntax-table)
  (run-hooks 'petscfem-mode-hook 'paragraph-indent-petscfem-mode-hook))

(defun paragraph-indent-minor-mode ()
  "Minor mode for editing text, with leading spaces starting a paragraph.
In this mode, you do not need blank lines between paragraphs when the
first line of the following paragraph starts with whitespace, as with
`paragraph-indent-mode'.
Turning on Paragraph-Indent minor mode runs the normal hook
`paragraph-indent-petscfem-mode-hook'."
  (interactive)
  (set (make-local-variable 'paragraph-start)
       (default-value 'paragraph-start))
  (set (make-local-variable 'paragraph-separate) paragraph-start)
  (run-hooks 'paragraph-indent-petscfem-mode-hook))
      
(defalias 'indented-petscfem-mode 'petscfem-mode)

(defun petscfem-mode-hook-identify ()
  "Mark that this mode has run `petscfem-mode-hook'.
This is how `toggle-petscfem-mode-auto-fill' knows which buffers to operate on."
  (make-local-variable 'petscfem-mode-variant)
  (setq petscfem-mode-variant t))

(add-hook 'petscfem-mode-hook 'petscfem-mode-hook-identify)

(defun toggle-petscfem-mode-auto-fill ()
  "Toggle whether to use Auto Fill in Text mode and related modes.
This command affects all buffers that use modes related to Text mode,
both existing buffers and buffers that you subsequently create."
  (interactive)
  (let ((enable-mode (not (memq 'turn-on-auto-fill petscfem-mode-hook)))
	(buffers (buffer-list)))
    (if enable-mode
	(add-hook 'petscfem-mode-hook 'turn-on-auto-fill)
      (remove-hook 'petscfem-mode-hook 'turn-on-auto-fill))
    (while buffers
      (with-current-buffer (car buffers)
	(if petscfem-mode-variant
	    (auto-fill-mode (if enable-mode 1 0))))
      (setq buffers (cdr buffers)))
    (message "Auto Fill %s in Text modes"
	     (if enable-mode "enabled" "disabled"))))

(defun center-paragraph ()
  "Center each nonblank line in the paragraph at or after point.
See `center-line' for more info."
  (interactive)
  (save-excursion
    (forward-paragraph)
    (or (bolp) (newline 1))
    (let ((end (point)))
      (backward-paragraph)
      (center-region (point) end))))

(defun center-region (from to)
  "Center each nonblank line starting in the region.
See `center-line' for more info."
  (interactive "r")
  (if (> from to)
      (let ((tem to))
	(setq to from from tem)))
  (save-excursion
    (save-restriction
      (narrow-to-region from to)
      (goto-char from)
      (while (not (eobp))
	(or (save-excursion (skip-chars-forward " \t") (eolp))
	    (center-line))
	(forward-line 1)))))

(defun center-line (&optional nlines)
  "Center the line point is on, within the width specified by `fill-column'.
This means adjusting the indentation so that it equals
the distance between the end of the text and `fill-column'.
The argument NLINES says how many lines to center."
  (interactive "P")
  (if nlines (setq nlines (prefix-numeric-value nlines)))
  (while (not (eq nlines 0))
    (save-excursion
      (let ((lm (current-left-margin))
	    line-length)
	(beginning-of-line)
	(delete-horizontal-space)
	(end-of-line)
	(delete-horizontal-space)
	(setq line-length (current-column))
	(if (> (- fill-column lm line-length) 0)
	    (indent-line-to
	     (+ lm (/ (- fill-column lm line-length) 2))))))
    (cond ((null nlines)
	   (setq nlines 0))
	  ((> nlines 0)
	   (setq nlines (1- nlines))
	   (forward-line 1))
	  ((< nlines 0)
	   (setq nlines (1+ nlines))
	   (forward-line -1)))))

;;; petscfem-mode.el ends here
