#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Ly=1;				# Long in transverse direction
$yratio = 10;			# refinement along y
$Ny = 20;			# Nbr of points along y
$Nx = 2*$Ny;			# Nbr of points along x
$Lx = 4;			# Long of comp. domain along x.
$Lslip=1;			# Entry slip lenght
$Lplate=1;			# Plate length
$abso = 0;			# Use absorbing b.c.'s?
$Twall = 0;			# Use fixed Twall condition?

$Machin = 0.5;			# Inlet Mach number

$gamma = 1.4;			# Cp/Cv ratio
$Rgas = 1;			# Gas constant number
$rhoref = 1;			# non-perturbed density 
$Tref = 1;			# wall temperature and reference value

@vars = qw(Rgas Nx Ny Lx Ly yratio Lslip Lplate
	   Machin gamma Rgas rhoref Tref
	   abso Twall);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkplate.m";

$n1 = ($Ny+1)*$Nx+1;
$n2 = ($Ny+1)*($Nx-1)+1;

1;
