#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

## mass-line geometry and position
$Ly = 1;			# length in `y' direction
$N = 40;			# Nbr of elements in `y' direction

$hy = $Ly/$N;			# mesh size in `y' direction
$Lx = $Ly;			# make square elements

$U = 2;
$tend = 12;
$Courant = 1;
$Dt = $Courant*$hy/$U;

$nstep = ceil($tend/$Dt);
$nsaverot = 1;
$nfile = -1;

# $adv_case = 'gaussian';
$adv_case = 'gaussian_diag';

@vars = qw(Ly N Lx hy  Courant U tend 
	   Dt nstep nsaverot nfile T adv_case); 
octave_export_vars(">data.m.tmp",@vars);
transcript2(@vars);             # print variables on output

if ($mkmesh==1) {
    system "octave -qH mkmesh.m > mkmesh.log.tmp";
}

1;
