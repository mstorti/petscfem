#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$Lx = 4;
$Nx = 100;

$perx = 0;			# solution is periodic 
$pery = 1;			# U(i+pery,j+perx) = U(i,j) 

$Machin = 0.5;
$gamma = 1.4;
$Rgas = 1;
$rhoref = 1;
$Tref = 1;

@vars = qw(Lx Nx perx pery Machin
	   gamma Rgas rhoref Tref);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkgfmovshock.m";
system "echo -n > gfmovshock.some-rslt.tmp";

1;
