#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 

$L = 1;				# Side of the cavity
$N = 20;
$hratio = 10 unless $hratio;			# refinement (h_center/h_wall)
$Re = 1000 unless $Re;		# Reynolds number
$viscosity = $L/$Re;		# viscosity
$Dt = 0.01;			# time step
$alehook = "$ENV{'ALEHOOK'}";
$uini = 0.2;

@vars = qw(L N Re viscosity Dt uini);
octave_export_vars(">data.m.tmp",@vars);

system("octave -qH mksqcav.m");

1;
