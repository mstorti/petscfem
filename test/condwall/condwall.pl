# $Id: condwall.pl,v 1.1 2005/03/28 02:27:45 mstorti Exp $ 

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl"; # Initializes ePerl 

$Lx1=1;
$Lx2=3;
$Ly=1;

$Nx1 = 10;
$Nx2 = 30;
$Ny = 10;
$hratio1 = 5;
$hratio2 = 10;

$DP = 1;

@vars = qw(Lx1 Lx2 Ly Nx1 Nx2 Ny hratio1 hratio2 DP);

octave_export_vars(">data.m.tmp",@vars);
system "octave -qH mkcondwall.m";

1;
