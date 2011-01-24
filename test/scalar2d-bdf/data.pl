#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

## mass-line geometry and position
$Ly = 1;			# length in `y' direction
$N = 20;			# Nbr of elements in `y' direction

$hy = $Ly/$N;			# mesh size in `y' direction
$Lx = $Ly;			# make square elements

$U = 2;
$tend = 1;
$Courant = 1;
$Dt = $Courant*$hy/$U;
$diffusive_jacobians = 0;

$nstep = ceil($tend/$Dt);
$nsaverot = 1;
$nfile = -1;

if ($adv_case eq 'bdf_dgcl') {
    $diffusive_jacobians = 1;
}

@vars = qw(Ly N Lx hy  Courant U tend 
	   Dt nstep nsaverot nfile T adv_case); 
octave_export_vars(">data.m.tmp",@vars);
transcript2(@vars);             # print variables on output

1;
