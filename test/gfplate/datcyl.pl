#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$R = 1;
$Rext = 5;

$Nr = 40;
$Nphi = 2*$Nr;
$rratio = 5;

$Machin = 0.7;
$gamma = 1.4;
$Rgas = 1;
$rhoref = 1;
$Tref = 1;
$abso = 1;
$restart = 0;
$use_symm = 0;

$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;

@vars = qw(R Rext Nr Nphi rratio Machin 
	   gamma Rgas rhoref Tref 
	   abso pref cref uref use_symm);
octave_export_vars(">data.m.tmp",@vars);

$oscript = ($use_symm ? 'mkcyl' : 'mkcyl2');
system "octave -qH $oscript.m > $oscript.log.tmp";

if (!$dx && !$restart) { system "echo -n > cylabso.some-rslt.tmp"; }

1;
