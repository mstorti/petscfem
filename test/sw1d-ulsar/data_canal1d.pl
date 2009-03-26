#! /usr/bin/perl
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
$Lx=2*10;
$Nx = 400;
$ulsar = 1;
$uref = 0;
$ndof=2;

if ($uref==1) {$nodele = 3;}
else {$nodele = 2;}

$Courant_fac = 1; 		# Dt is computed from this
$Dt = 0.05*10;
$restart = 0;
$nstep = 2000;
$gravity = 1.;
#initial conditions
$u_ini = 0.5;
$h_ini = 1.0;

#inlet conditions
$u_in = $u_ini;
$h_in = $h_ini;

#outlet conditions
$u_out = $u_ini;
$h_out = $h_ini;

$Uref=[$u_out, $h_out];
## related to the vertical
$wall_angle = 45; ## triang y trape

$angle_ap = 270; ## circular2

$B1 = 4;
$B2 = 6;
$Z1 = 2;

## bump data
$mu = $Lx/2;
$sig = 0.04*$Lx;
$rc = 0.1*$Lx;
$A = 2.;

@vars = qw(ndof Courant_fac Dt restart nstep gravity  
          h_ini u_ini u_in h_in u_out h_out
          Nx Lx ulsar Uref wall_angle uref nodele
	  mu sig rc A B1 B2 Z1);

octave_export_vars(">data.m.tmp",@vars);

1;
