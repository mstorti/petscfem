
run:
	octave -qH mksqcav.m
	mpirun -np 1							     \
		/u/mstorti/PETSC/petscfem-beta-3.30/applications/ns/ns_O.bin \
		-case sqcav.depl
