#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 50;

$phi = 0;
$theta = 0;

$Machin = 0.5;
$gamma = 1.4;
$Rgas = 1;
$rhoref = 1;
$Tref = 1;
$dufac = 0.01;
$sigma = 0.3;
$abso = 1;

$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;
$Uref = [$rhoref,$uref,0,0,$pref];
$norx = cos($theta)*cos($phi);
$nory = cos($theta)*sin($phi);
$norz = sin($theta);
$nor = [$norx,$nory,$norz];

@vars = qw(sigma Rgas Nx Lx Machin gamma
	   Rgas rhoref Tref dufac abso uref
	   pref cref Uref phi theta norx nory 
	   norz nor longindx);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso.m";

1;
