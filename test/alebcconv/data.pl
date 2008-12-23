#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 1;
$Nx = 50;
$Courant = 0.5;

## W.r.t. the walls 
$Machin = 0.5;
$ga = 1.4;
$rhoref = 1;
$pref = 1/$ga;
$use_ALE = 1 unless defined $use_ALE;

$cref = sqrt($ga*$pref/$rhoref);
$uref = $Machin*$cref;
$Uref = [$rhoref,$uref,0,$pref];

$Dt = $Lx/$Nx/($cref+$uref);

$ufl = $uref;
$vmesh = ($use_ALE ? -$uref : 0);
$ufl += $vmesh;

@vars = qw(Lx Nx Dt Machin ga Uref cref 
	   uref rhoref pref use_ALE use_linadv use_invbc
	   ufl vmesh use_perix du);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso.m";

1;
