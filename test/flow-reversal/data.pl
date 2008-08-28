#$Id: strip2d.epl,v 1.4 2003/01/10 16:28:52 mstorti Exp $
# Initializes ePerl 
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
#
$use_exterior_normal = 1;
$Lx = 2;
$Ly = 1;
$Nx=$Ny=10;
$gb = 1.0 unless defined $gb;
@vars = qw(Lx Ly Nx Ny use_exterior_normal gb);
transcript("", @vars);	# print variables on output and transcript this block
octave_export_vars(">data.m.tmp",@vars);
system "octave -qH mkmesh.m" if $mkmesh;

1;
