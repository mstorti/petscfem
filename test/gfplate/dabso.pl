#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx=4;
$Nx = 80;

$Machin = 0.5;
$gamma = 1.4;
$Rgas = 287;
$rhoref = 1;
$Tref = 300;

@vars = qw(Rgas Nx Lx Machin gamma Rgas rhoref Tref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso.m";

1;
