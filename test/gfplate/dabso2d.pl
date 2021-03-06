#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 50;
$rota = 1;

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
$Uref = [$rhoref,$uref,0,$pref];
$alpha = $PI/2*$rota;		# Rotation angle
$longindx = 1+($rota % 2);	# Longitudinal indx
$norx = cos($alpha);
$nory = sin($alpha);

@vars = qw(sigma Rgas Nx Lx Machin gamma
	   Rgas rhoref Tref dufac abso uref
	   pref cref Uref alpha longindx);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso2d.m";
system "echo -n > gfabso2d.some-rslt.tmp";

1;
