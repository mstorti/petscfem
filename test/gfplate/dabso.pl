#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 50;

$Machin = -0.5;
$gamma = 1.4;
$Rgas = 287;
$rhoref = 1;
$Tref = 300;
$dufac = 0.01;
$sigma = 0.3;
$abso = 1;

$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;
$Uref = [$rhoref,$uref,0,$pref];

@vars = qw(sigma Rgas Nx Lx Machin gamma
	   Rgas rhoref Tref dufac abso uref
	   pref cref Uref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso.m";

1;
