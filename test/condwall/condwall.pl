# $Id: condwall.pl,v 1.5 2005/03/29 12:47:01 mstorti Exp $ 

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl"; # Initializes ePerl 

$Lx1=1;
$Lx2=1;
$Ly=1;

$Nx1 = $Nx2 = $Ny = 40;
$hratio1 = 5;
$hratio2 = 5;

$DP = 1;
$use_peri = 0;
$Lslit = $Ly/2;

@vars = qw(Lx1 Lx2 Ly Nx1 Nx2 Ny hratio1 hratio2 DP
	   use_peri Lslit);

octave_export_vars(">data.m.tmp",@vars);
system "octave -qH mkcondwall.m";

1;
