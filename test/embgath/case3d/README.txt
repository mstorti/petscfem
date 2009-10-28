The test is a 3D strip of N elements occupying 
0<= x,z <= h, 0<=y<=L, with h=L/N, and N is the number of elements in
the y direction. 

The elements are hexas. 

Periodic boundary conditions are imposed in x and z directions. 

A body force G_body drives the flux in the x direction. 

The solution is u(x) = Umax *(4/L^2)*x*(L-x), so that the momentum
balance is

rho*G_body = (8/L^2)*mu*Umax

We make various checks changing the values of rho, visco(=mu), L,
G_body. 

The torque is measured with respect to the (h/2,0,h/2). 

WARNING: this tests are somewhat coincident with tests in
test/nutatcyl 
