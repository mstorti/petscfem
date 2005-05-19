PETSCFEM_DIR := /u/mstorti/PETSC/petscfem-beta-3.30

run:
	octave -qH mksqcav.m
	mpirun -np 1							     \
		$(PETSCFEM_DIR)/applications/ns/ns_O.bin \
		-case sqcav.depl
