#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 100;
$rota = 0;

$Machin = 0.0;
$gamma = 1.4;
$Rgas = 1;
$rhoref = 1;
$Tref = 1;
$dufac = 0.01;
$sigma = 0.3;
$abso = 0;
$all_fields = 1;

$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;
$Uref = [$rhoref,$uref,0,$pref];
$alpha = $PI/2*$rota;		# Rotation angle
$norx = cos($alpha);
$nory = sin($alpha);

$longindx = $rota % 2 + 1;

@vars = qw(sigma Rgas Nx Lx Machin gamma
	   Rgas rhoref Tref dufac abso uref
	   pref cref Uref alpha longindx
	   all_fields);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso2dn.m";
system "echo -n > gfabso2dn.some-rslt.tmp";

1;
