# -*- mode: shell-script -*-
#
<:require '../eperlini.pl':>// # Initializes ePerl 

# Synchronize the time step with external frequency
<:
$Nx = 3;
$Ny = 2;
@vars = qw(Nx Ny);
octave_export_vars(">data.m.tmp",@vars);
:>//

global_options

nsave 100000000
save_file oscplate.sal
nstep 5

#    Iterative parameters
nnwt 3
tol_newton 0
atol 0
rtol 0
dtol 1e3
solver "petsc"
preco_type "lu"

initial_state "couette.ini.tmp"
Dt 1
steady
alpha 1.
G_body 1. 0.
tau_fac 1.
print_linear_system_and_stop
solve_system 0
verify_jacobian_with_numerical_one
# inwt_stop 2
__END_HASH__

# ndim nu ndof
nodes  2  2  8
__INCLUDE__ couette.nod.tmp
_END_NODES__

elemset bubbly 4
props 
#
#Datos del elemento
#
name malla
geometry cartesian2d
ndim 2
npg 4
#
# Datos fisicos
#
weak_form 1
C_mu 0.
rho_l 1.
rho_g 1e-3
visco_l 1.1e-3
visco_g 1e-3
__END_HASH__
__INCLUDE__ couette.con.tmp
__END_ELEMSET__

end_elemsets

fixa

#if 0
# Pressure at some node
1 1 0.

# velocity set to zero at the bottom
1 3 0.
1 4 1.
2 3 0.
2 4 1.
<:=2*$N+1:> 3 0.
<:=2*$N+1:> 4 0.
<:=2*$N+2:> 3 0.
<:=2*$N+2:> 4 0.
#else
#if 0
1 1 0.
1 3 0.
1 4 0.
2 1 0.
2 3 0.
2 4 0.
3 1 0.
3 3 0.
3 4 0.
4 1 0.
4 3 0.
4 4 0.
6 1 0.
6 3 0.
6 4 0.
7 1 0.
7 3 0.
7 4 0.
9 1 0.
9 3 0.
9 4 0.
10 1 0.
10 3 0.
10 4 0.
11 1 0.
11 3 0.
11 4 0.
12 1 0.
12 3 0.
12 4 0.
#endif
#endif
__INCLUDE__ couette.fixa.tmp
__END_FIXA__

#if 0
#
# Periodic boundary conditions from y=0.1 to y=0.
#
constraint
__INCLUDE__ couette.peri.tmp
__END_CONSTRAINT__
#endif
