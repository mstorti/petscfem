#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 1.797525440;
$Nx = 10;
$Nr = 5;
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

if (0) { 
    $pout = 143; 
    $Tout = 262;
} else { 
## 0.00023833 is the target $pin/$pout
#    $pout = 143;
    $pout = 0.7*$pin;
#    $Tout = 262;
    $Tout = $Tin;
    $rhoout = $pout/($Rgas*$Tout);
}
$pout0 = $pout/$pref;
$rhoout0 = $rhoout/$rhoref;

$Co = 0.5;
$h = 1/$Nx;
$Dt = $Co*$h/($uin0+1);
$tramp = 20*$Dt;

@vars = qw(Nr Dt ga Lx Nx Rgas pin0 rhoin0 uin0
	   pout0 rhoout0 tramp dtheta);
octave_export_vars(">data.m.tmp",@vars);
doc_vals(@vars);

$mkmesh = $ENV{'mkmesh'};
$mkmesh = 1 if !defined $mkmesh;
if (!defined $mkmesh || $mkmesh) {
    system "octave -qH mkgfshock3d.m";
}
system "echo -n > gfshock.some-rslt.tmp";

1;
