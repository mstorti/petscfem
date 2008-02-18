#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
$Nx = 100;
$Ny = 20;
$turbulent = 0;

if ($turbulent) {$ndof= 5;}
else {$ndof=3;}

$Courant_fac = 1; 		# Dt is computed from this
$Dt = 0.05;
$restart = 0;
$nstep = 200;
$gravity = 9.81;
#initial conditions
$u_ini = 0.5;
$v_ini = 0.0;
$h_ini = 1.0;

#inlet conditions
$u_in = 0.5;
$v_in = 0.0;
$h_in = 1.0;

#outlet conditions
$u_out = 0.0;
$v_out = 0.025;
$h_out = 0.95;

$kappa = 0.1;
$epsilon = 0.1;

@vars = qw(ndof Courant_fac Dt restart nstep gravity C_mu 
          h_ini u_ini v_ini u_in v_in h_in u_out v_out h_out
          kappa epsilon Nx Ny turbulent);

octave_export_vars(">data.m.tmp",@vars);

#system "octave -qH mkplanta.m";

1;

