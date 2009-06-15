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
elsif ($subcase eq 'rho10') { $rho *= 10.0; }
elsif ($subcase eq 'g10') { $gbody *= 10.0; }
elsif ($subcase eq 'l10') { $L *= 10.0; }
elsif ($subcase eq 'all10') { 
    $L *= 10.0; 
    $visco *=10.0;
    $gbody *=10.0;
    $rho *=10.0;
}
else { die "unknown case $case"; }

$h = $L/$N;
@vars = qw(L h hratio N visco rho gbody subcase);

# print variables on output and transcript this block
transcript("", @vars);

octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkmesh.m >mkmesh.out.tmp" if $mkmesh;

1;
