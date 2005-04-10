#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Ly=1;				# Long in transverse direction
$DLy=0.2;  			# Reduction of aperture at center
$yratio = 1;			# refinement along y
$Ny = 10;			# Nbr of points along y
$Nx1 = 2*$Ny;			# Nbr of points along x (duct)
$Nx2 = $Ny;			# Nbr of points along x (exterior)
$Lx1 = 2*$Ly;			# Long of comp. domain along x. (duct)
$Lx2 = 3*$Ly;			# Long of comp. domain along x. (exterior)
$Rgas = 1;
$Tref = 1;
$rhoref = 1;
$gamma = 1;
$Machin = 0.5;
$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;
$Uref = [$rhoref,$uref,0,$pref];

@vars = qw(Machin gamma rhoref Tref Rgas
	   Ly DLy yratio Ny Nx1 Nx2  Lx1 Lx2 
	   pref cref uref Uref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mknozzle2.m";

1;

