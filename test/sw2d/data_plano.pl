#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
$Lx=98;
$Ly=1;
$Nx = 1000;
$Ny = 10;
$turbulent = 0;

if ($turbulent) {$ndof= 5;}
else {$ndof=3;}

$Courant_fac = 1; 		# Dt is computed from this
$Dt = 0.05;
$restart = 0;
$nstep = 10000;
$gravity = 9.81;
#initial conditions
$B=$Ly;
$q=0.45;

$h_ini = 0.35;
$u_ini = 1.;#$q/($B*$h_ini);
$v_ini = 0.0;

#inlet conditions
$h_in = 0.46;
$u_in = $q/($B*$h_in);
$v_in = 0.0;

#outlet conditions
$h_out = 0.3;
$u_out = $q/($B*$h_out);;
$v_out = 0.0;

$kappa = 0.1;
$epsilon = 0.1;

@vars = qw(ndof Courant_fac Dt restart nstep gravity C_mu 
          h_ini u_ini v_ini u_in v_in h_in u_out v_out h_out
          kappa epsilon Lx Ly Nx Ny turbulent);

octave_export_vars(">data.m.tmp",@vars);

1;

