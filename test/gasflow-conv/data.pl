# $Id: cubcav.pl,v 1.8 2006/01/27 20:51:51 mstorti Exp $
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
require "math.pl";

$check_conv = 1 unless defined $check_conv;

$n = 1;
$L = 1;
$weak_form = 1;
$Machin = 2;
$Co = 1;
$peri_x = 1;

$nu = 1e-10;
$gamma = 1.4;
$pref = 1.0/$gamma;
$rhoref = 1.0;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;
## Sol for M=2,  U2= [2.66667   0.75000   3.21429]
$p2 = 3.21429;

$h = $L/$n;
$Dt = $h/($uref+$cref);

$T = 3.0*(0.5*$L)/$cref;
$nstep = int(ceil($T/$Dt));

@vars= qw(Dt n nu L weak_form gamma pref cref rhoref uref p2);
octave_export_vars(">data.m.tmp",@vars);
transcript2(@vars);	# print variables on output

system "octave -qH mkmesh.m > mkmesh.output.tmp" if $mkmesh;

1;
