(define my-rank (mpi-rank))
(define size (mpi-size))
(format #t "myrank ~A, size ~A\n" my-rank size)

(define val #f)

#!
(do ((j 0 (+ j 2))) ((= j 20)) 
  (cond ((= my-rank 0)
	 (mpi-send j 1)
	 (set! val (mpi-recv 1)))
	(#t
	 (set! val (mpi-recv 0))
	 (mpi-send (+ j 1) 0)))
  (format #t "[~A] received ~A\n" my-rank val))
!#

(do ((j 0 (+ j 1))) ((= j 20)) 
  (cond ((= my-rank 1)
	 (mpi-send 23 0))
	(#t
	 (set! val (mpi-recv 1))
	 (format #t "in guile: received ~A\n" val))))

(mpi-finalize)
