#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
$Lx=2;
$Ly=0.2;
$Nx = 400;
$Ny = 1;
$ulsar = 0;
$turbulent = 1;

if ($ulsar==1) {$nodele = 2;}
else {$nodele=3;}
 
if ($turbulent) {$ndof= 5;}
else {$ndof=3;}

$Courant_fac = 1; 		# Dt is computed from this
$Dt = 0.005;
$restart = 0;
$nstep = 1000;
$gravity = 9.81;
#initial conditions
$u_ini = 0.;
$v_ini = 0.;
$h_ini = 1.0;

#inlet conditions
$u_in = $u_ini;
$v_in = $v_ini;
$h_in = $h_ini;

#outlet conditions
$u_out = $u_ini;
$v_out = $v_ini;
$h_out = $h_ini;

$kappa = 0.1;
$epsilon = 0.1;

#$Uref=[$u_out, $v_out, $h_out, $kappa, $epsilon];

@vars = qw(ndof Courant_fac Dt restart nstep gravity C_mu 
          h_ini u_ini v_ini u_in v_in h_in u_out v_out h_out
          kappa epsilon Nx Ny turbulent Lx Ly ulsar nodele);

octave_export_vars(">data.m.tmp",@vars);

1;

