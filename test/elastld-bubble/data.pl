#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

## mass-line geometry and position
$ref = 2;
$R = 1;
$DR = 0.1;
$Nr = 2*$ref;
$Nphi = 20;
$Dt = 0.01;
$Tend = 100;
$Tstop = -1;
$nstep = ceil($Tend/$Dt);
$Eratio = 1;

@vars = qw(ref R  DR Nr Nphi Dt Tend Tstop nstep Eratio); 
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkmesh.m > mkmesh.log.tmp";

1;
