#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

## mass-line geometry and position
$ref = 1;
$Lx = 0.5;
$Nx = 2*$ref;
$Ly = 5;
$Ny = 10*$ref;
$ndim = 3;
$Rratio = 0.5; # in 3D: Rext = Lx/2, Rint= Rratio*Rext
$Nphi = 10*$ref;
$Dt = 0.05;
$Tend = 2.0;
$Tstop = -1;
$nstep = ceil($Tend/$Dt);
$Eratio = 1;
$ntrail = 100;

if ($ndim==2) {
    $nel = 4;
    $npg = 4;
} elsif ($ndim==3) {
    if ($use_tetra) {
        $nel = 4;
        $npg = 4;
    } else {
        $nel = 8;
        $npg = 8;
    }
} else {
    die "bad ndim $ndim";
}

@vars = qw(Lx Ly Nx Ny Nphi ndim nel npg Rratio Dt 
           Tend Tstop nstep Eratio Ntrail); 
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkmesh.m > mkmesh.log.tmp";

1;
