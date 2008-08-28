#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
$Lx=1;
$Ly=0.2;
$Nx = 200;
$Ny = 1;
$ulsar = 0;
$abso = 0;
if ($ulsar==1) {$nodele = 2;}
else {$nodele=3;}
 
$restart = 0;
$nstep = 1000;
$gravity = 1;
$rho1 = 1.;
$rho2 = 0.5;

#initial conditions
$u1_ini = 0.4;
$v1_ini = 0.;
$h1_ini = 1.0;
$u2_ini = 0.5;
$v2_ini = 0.;
$h2_ini = 1.0;

#inlet conditions
$u1_in = $u1_ini;
$v1_in = $v1_ini;
$h1_in = $h1_ini;
$u2_in = $u2_ini;
$v2_in = $v2_ini;
$h2_in = $h2_ini;

# $u1_in = $u1_ini;
# $v1_in = $v1_ini;
# $h1_in = $h1_ini;
# $u2_in = $u2_ini;
# $v2_in = $v2_ini;
# $h2_in = $h2_ini;

#outlet conditions
$u1_out = $u1_ini;
$v1_out = $v1_ini;
$h1_out = $h1_ini;
$u2_out = $u2_ini;
$v2_out = $v2_ini;
$h2_out = $h2_ini;

$ndof = 6;
$hx = $Lx/$Nx;
$Courant_fac = 1; 		# Dt is computed from this
$Dt = $hx*$Courant_fac/$u2_ini;

@vars = qw(ndof Courant_fac Dt restart nstep gravity rho1 rho2
          h1_ini u1_ini v1_ini u1_in v1_in h1_in u1_out v1_out h1_out
          h2_ini u2_ini v2_ini u2_in v2_in h2_in u2_out v2_out h2_out
          Nx Ny Lx Ly ulsar nodele abso hx);

octave_export_vars(">data.m.tmp",@vars);

1;

