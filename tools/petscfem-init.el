;;; petscfem-init.el --- Initialize PETSc-FEM mode

;; Copyright (C) 2003  Free Software Foundation, Inc.

;; Author:  <mstorti@spider>
;; Keywords: 

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

;; These are common tasks you probably want to put in your .emacs
;; so that PETSc-FEM mode is correctly initialized. 

;;; Code:

;;; $Id: petscfem-init.el,v 1.2 2003/11/24 00:06:50 mstorti Exp $

;; Load info-look if not already loaded. 
(load-library "info-look")

;; Add PETSc-FEM mode to the list of modes supported by info-look
;; You can replace this by a straight path, i.e. 
;; (defvar petscfem-info-file "/home/bob/petscfem/options.info")

;; Configure this if needed. Not necessary if PETSCFEM_DIR
;; environment variable already defined. 
; (setq petscfem-dir "/home/bob/petscfem")

;; Configure this if needed. Not necessary if `petscfem-dir'
;; is defined and want to point to the  `tools/options.info'
; (setq petscfem-info-file "/home/bob/petscfem/tools/options.info")

;; Load `info-look' definitions for `info-lookup'
(info-lookup-maybe-add-help
 :mode 'petscfem-mode 
 :regexp "[_a-zA-Z0-9./+-]+"
 :doc-spec (list (list (concat petscfem-info-file "Options Index")))

;; This key bindings are util for pasting option names into the
;; data file. See the Emacs section in the PETSc-FEM documentation. 
(define-key Info-mode-map (kbd "c") 'my-Info-lookup-copy-keyword)
(define-key Info-mode-map (kbd "x") 'my-Info-bury-and-kill)

(provide 'petscfem-init)
;;; petscfem-init.el ends here
