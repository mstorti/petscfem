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

;;; $Id: petscfem-init.el,v 1.1 2003/11/23 23:36:33 mstorti Exp $

;;; Commentary:

;; These are common tasks you probably want to put in your .emacs
;; so that PETSc-FEM mode is correctly initialized. 

;;; Code:

;; Load info-look if not already loaded. 
(load-library "info-look")
;; Add PETSc-FEM mode to the list of modes supported by info-look
(info-lookup-maybe-add-help
 :mode 'petscfem-mode 
 :regexp "[_a-zA-Z0-9./+-]+"
 :doc-spec (list (list 
		  (concat "(" (getenv "PETSCFEM_DIR") 
			  "/doc/options.info)Options Index"))))
(define-key Info-mode-map (kbd "c") 'my-Info-lookup-copy-keyword)
(define-key Info-mode-map (kbd "x") 'my-Info-bury-and-kill)

(provide 'petscfem-init)
;;; petscfem-init.el ends here
