#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Rgas = 287;
$Nx=40;
$Ny=20;
$Lx = 4;
$Lplate=1;

$Machin = 0.5;

$gamma = 1.4;
$Rgas = 287;
$rhoref = 1;
$Tref = 300;

@vars = qw(Rgas Nx Ny Lx Lplate Machin gamma Rgas rhoref Tref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkplate.m";

1;
