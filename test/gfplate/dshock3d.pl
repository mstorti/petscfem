#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 1.797525440;
$Nx = 100;
$Nr = 1;
$dtheta = 1*$PI/180;

$ga = 1.17;
$Molw = 27.5;
$Molw_air = 29;
$Rgas_air = 287;
$Rgas = $Rgas_air*$Molw_air/$Molw;

$pin = 6e5;
$Tin = 4170;
$rhoin = $pin/($Rgas*$Tin);
$uin = 0;

$uref = sqrt($ga*$pin/$rhoin);
$rhoref = $rhoin;
$pref = $rhoref*$uref**2.;

$pin0 = $pin/$pref;
$rhoin0 = $rhoin/$rhoref;
$uin0 = $uin/$uref;

$pout = 0.01*$pin;
$Tout = 262;
$rhoout = $pout/($Rgas*$Tout);

$pout0 = $pout/$pref;
$rhoout0 = $rhoout/$rhoref;

$Co = 0.5;
$h = $Lx/$Nx;
$Dt = $Co*$h/($uin0+1);
$tramp = 20*$Dt;

$Rscale = 0.02;

@vars = qw(Nr Dt ga Lx Nx Rgas pin0 rhoin0 uin0
	   pout0 rhoout0 tramp dtheta Rscale);
octave_export_vars(">data.m.tmp",@vars);
doc_vals(@vars);

$mkmesh = $ENV{'mkmesh'};
$mkmesh = 0 if !defined $mkmesh;
if (!defined $mkmesh || $mkmesh) {
    system "octave -qH mkgfshock3d.m";
}

system "echo -n > gfshock3d.some-rslt.tmp"
    unless $dx;

1;
