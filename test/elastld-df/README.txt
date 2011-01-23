This test checks the large displacement elasticity elemset
ld_elasticity. Previously we had a formulation based on displacements
and velocities, in order to have a first order PDE in time, so as to
use the trapezoidal rule. But now we implemented the Newmark method in
a sepcial main `struct_main()' (file nsstruct.cpp). 

This test is a hollow cylinder where the bottom experiments a
sinusoidal movement, with a frequency close to the first natural mode
of the structure. We do this for 2 seconds, whereas the period of
excitation is 2.5 seconds. We compare the result with a previously
computed reference solution in elastld.state40.ref. 
