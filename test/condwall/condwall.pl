# $Id: condwall.pl,v 1.2 2005/03/28 16:42:59 mstorti Exp $ 

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl"; # Initializes ePerl 

$Lx1=1;
$Lx2=1;
$Ly=1;

$Nx1 = 20;
$Nx2 = 20;
$Ny = 1;
# $hratio1 = 5;
# $hratio2 = 10;
$hratio1 = 1;
$hratio2 = 1;

$DP = 1;
$use_peri = 0;

@vars = qw(Lx1 Lx2 Ly Nx1 Nx2 Ny hratio1 hratio2 DP
	   use_peri);

octave_export_vars(">data.m.tmp",@vars);
system "octave -qH mkcondwall.m";

1;
