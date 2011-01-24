This test checks the large displacement elasticity elemset
ld_elasticity. Previously we had a formulation based on displacements
and velocities, in order to have a first order PDE in time, so as to
use the trapezoidal rule. But now we implemented the Newmark method in
a sepcial main `struct_main()' (file nsstruct.cpp). 

This test is an annulus  |r-R|<DR, with R=1, DR=0.1. The annulus is
imposed an external differential pressure of 0.1 By symmetry only a
sector 0 <= \phi <= \pi/2 is considered. 
