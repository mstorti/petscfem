#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 80;

$Machin = 0.5;
$gamma = 1.4;
$Rgas = 287;
$rhoref = 1;
$Tref = 300;
$du = 0.3;
$sigma = 0.3;

@vars = qw(sigma Rgas Nx Lx Machin gamma Rgas rhoref Tref du);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfabso.m";

1;
