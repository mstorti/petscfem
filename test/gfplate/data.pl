#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Ly=1;
$yratio = 10;
$ref = 1;
$Ny = 15;
$Nx = 2*$Ny;
$Lx = 4;
$Lplate=1;

$Machin = 0.5;

$gamma = 1.4;
$Rgas = 1;
$rhoref = 1;
$Tref = 1;

@vars = qw(Rgas Nx Ny Lx Ly yratio Lplate Machin gamma Rgas rhoref Tref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkplate.m";

$n1 = ($Ny+1)*$Nx+1;
$n2 = ($Ny+1)*($Nx-1)+1;

1;
