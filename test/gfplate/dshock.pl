#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 50;

$gamma = 1.17;
$Molw = 27.5;
$Molw_air = 29;
$Rgas_air = 287;
$Rgas = $Rgas_air*$Molw_air/$Molw;

$pin = 6e5;
$Tin = 4170;
$rhoin = $pin/$Rgas*$Tin;
$uin = 0;

$fac = 1000;
$pout = 143;
$Tout = 4170;
$rhoout = $pout/$Rgas*$Tout;
$uout = 0;

@vars = qw(Lx Nx Rgas pin rhoin uin pout rhoout uout);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfshock.m";
# system "echo -n > gfshock.some-rslt.tmp";

1;
