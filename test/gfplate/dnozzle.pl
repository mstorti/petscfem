#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Ly=1;				# Long in transverse direction
$yratio = 10;			# refinement along y
$Ny = 20;			# Nbr of points along y
$Nx = 2*$Ny;			# Nbr of points along x
$Lx = 4;			# Long of comp. domain along x.
$Aratio = 2; 			# Area at inlet/area ate outlet
$Lnozzle = 2;			# Contraction is done in this length
$Rgas = 1;
$Tref = 1;
$rhoref = 1;
$gamma = 1;
$Machin = 0.5;

@vars = qw(Machin gamma rhoref Tref Rgas
	   Ly yratio Ny Nx  Lx Aratio Lnozzle);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mknozzle.m";

1;

