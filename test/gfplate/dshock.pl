#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 1.8;
$Nx = 50;

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

if (0) { 
    $pout = 143; 
    $Tout = 262;
} else { 
    $pout = 0.99*$pin; 
    $Tout = 262;
    $rhoout = $pout/($Rgas*$Tout);
}
$pout0 = $pout/$pref;
$rhoout0 = $rhoout/$rhoref;

$Co = 0.2;
$h = 1/$Nx;
$Dt = $Co*$h/2;
$tramp = $100*$Dt;

@vars = qw(Dt ga Lx Nx Rgas pin0 rhoin0 uin0 pout0 rhoout0 tramp);
octave_export_vars(">data.m.tmp",@vars);
doc_vals(@vars);

$mkmesh = $ENV{'mkmesh'};
if (!defined $mkmesh || $mkmesh) {
    system "octave -qH mkgfshock.m";
}
system "echo -n > gfshock.some-rslt.tmp";

1;