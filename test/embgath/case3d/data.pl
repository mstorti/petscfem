#$Id: strip2d.epl,v 1.4 2003/01/10 16:28:52 mstorti Exp $
# Initializes ePerl 
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";
#
$N = 20;
$L = 1;
$hratio = 5;
@vars = qw(L hratio N);

# print variables on output and transcript this block

octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkmesh.m >mkmesh.out.tmp" if $mkmesh;

1;
