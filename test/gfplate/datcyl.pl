#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$R = 1;
$Rext = 5;

$Nr = 20;
$Nphi = $Nr;
$rratio = 5;

$Machin = 0.2;
$gamma = 1.4;
$Rgas = 1;
$rhoref = 1;
$Tref = 1;
$abso = 1;
$restart = 0;

$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;

@vars = qw(R Rext Nr Nphi rratio Machin 
	   gamma Rgas rhoref Tref 
	   abso pref cref uref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkcyl.m";
if (!$dx && !$restart) { system "echo -n > cylabso.some-rslt.tmp"; }

1;
