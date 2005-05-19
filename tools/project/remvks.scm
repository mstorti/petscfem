;;; ttn/optargs-kw-utils.scm --- make using (ice-9 optargs-kw) easier

;; Rel:v-0-41-pianto-due-ore
;;
;; Copyright (C) 2003-2005 Thien-Thi Nguyen
;; This file is part of ttn's personal scheme library, released under GNU
;; GPL with ABSOLUTELY NO WARRANTY.  See the file COPYING for details.

;;; Commentary:

;; This module exports:
;;  (remove-keys ls)

;;; Code:

(define-module (ttn optargs-kw-utils)
  #:export (remove-keys))

(define (remove-keys ls)
  (let loop ((ls ls) (acc '()))
    (if (null? ls)
        (reverse! acc)
        (let ((kw? (keyword? (car ls))))
          (loop ((if kw? cddr cdr) ls)
                (if kw? acc (cons (car ls) acc)))))))

;;; ttn/optargs-kw-utils.scm ends here
