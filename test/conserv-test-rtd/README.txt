This example was taken from the calculation of the Residence Time
Distribution (RTD) for particles. We have a container of shape `Omega'
with an Sinlet and Soutlet surfaces. There is a solenoidal velocity
vector, constant in time, and we inject a particles at a constant
density at the inlet. The RTD is the probability E(t) that a particle
exits at a time t. This can be computed by probabilistic,
i.e. MonteCarlo methods: we inject particles at the inlet, wait that
they exit and take note of this time. Then we compute an
histogram. Another possibility is to use the `tracer' method, i.e. we
consider the concentration of a tracer c. We assume thath initially
c=0 in Omega and then at t=0 we impose c=1 at Sinlet. Then we solve
for the advection diffusion PDE and compute the averaged value of c at
Soutlet. 

cout(t) = (\int_Soutlet c(x,t) vn d\Sigma)
           / (\int_Soutlet vn d\Sigma)

Note that the average must be weighted to the local normal
velocity. Then it can be shown by a balance of concentration in Omega
that

tmean = \int_0^\infty (1-cout(t)) dt

In this test we take a rectangular 2D duct of length Lz, and
`diameter' Ly=2*radius. The flow is assumed Poiseuille 

ux(y) = (1-(y/radius)^2)*Vmax

con Vmax = 1.5*Vmean. 

It can be shown that tmean= Lz/Vmean. So we set the velocity ux(y) to
the analytic value and solve the PDE for c(x,t) and then compute the
`tmean' and compare with the analytic value. 
