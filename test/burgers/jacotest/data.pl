require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
require "math.pl";

$check_jac = 1 unless defined $check_jac;
$N = ($check_jac? 1 : 50);
$steady = ($check_jac? 0 : 1);
$diffusivity = ($check_jac? 0.01 : 0.1);
$L=1;

@vars = qw(N L);
octave_export_vars(">data.m.tmp",@vars);
system "octave -qH mkmesh.m";

1;
