#$Id: strip2d.epl,v 1.4 2003/01/10 16:28:52 mstorti Exp $
# Initializes ePerl 
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
#
$N = 20;
$L = 1;
$hratio = 5;

$visco = 1.0;
$rho = 1.0;
$gbody = 1.0;

if ($subcase eq 'ref') {}
elsif ($subcase eq 'visco10') { $visco *= 10.0; }
else { die "unknown case $case"; }

@vars = qw(L hratio N visco rho gbody subcase);

# print variables on output and transcript this block
transcript("", @vars);

octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkmesh.m >mkmesh.out.tmp" if $mkmesh;

1;
