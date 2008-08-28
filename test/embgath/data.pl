#$Id: strip2d.epl,v 1.4 2003/01/10 16:28:52 mstorti Exp $
# Initializes ePerl 
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
#
$casename = $ENV{'CASE_NAME'};       # This is used in the `proc?.m' scripts
$ini = $ENV{'ini'};             # Run the initialization step
$use_exterior_normal = 1;
$gather = 0;
$props = 0;
$dump = 1;
$L = 1;
$Nx=$Ny=10;
$h=$L/$Nx;
$ndim = 2;
$Dt=0.1;          # Time step
$maxits = 100;
$xratio = 5;
$layers = 2;
@vars = qw(L Nx Ny h xratio use_exterior_normal);
transcript("", @vars);	# print variables on output and transcript this block
octave_export_vars(">data.m.tmp",@vars);
system "octave -qH mkstrip2d.m >mkstrip2d.out.tmp" if $mkmesh;

1;
