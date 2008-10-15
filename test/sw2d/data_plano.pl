#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
$L = 1.0;
$Nx = 100;
$Ny = 100;
$turbulent = 0;

if ($turbulent) {$ndof= 5;}
else {$ndof=3;}

$Courant_fac = 1;
$Dt = 0.005;
$restart = 0;
$nstep = 2000;
$gravity = 9.81;
#initial conditions
$u_ini = 0.0;
$v_ini = 0.0;
$h_ini = 1.0;

$A=10;
$B=$A;
$sigma=0.7;
$x0=$L/4;
$y0=$x0*1.5;
$ampli = 6.5;
$base = 1.0;
## 2D Gauss curve
## z=1+0.5.*(1./(2*pi).*exp(-a*((xx-0.5)./sigma).^2-b*((yy-0.75)./sigma).^2));

#inlet conditions
$u_in = 1.0;
$v_in = 0.0;
$h_in = 1.0;

#outlet conditions
$u_out = 0.0;
$v_out = 0.0;
$h_out = 1.0;

$kappa = 0.1;
$epsilon = 0.1;

@vars = qw(ndof Courant_fac Dt restart nstep gravity C_mu 
          h_ini u_ini v_ini u_in v_in h_in u_out v_out h_out
          kappa epsilon L Nx Ny turbulent
          A B sigma x0 y0 ampli base);

octave_export_vars(">data.m.tmp",@vars);

1;

