# $Id: cubcav.pl,v 1.4 2003/11/27 00:07:32 mstorti Exp $
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
require "math.pl";

$case = 'laplace' unless $case;
get_var_env2('data_dir','.');
$endur_test = $data_dir ne '.';
get_var_env2('job',undef);
get_var_env2('mkmesh',1);
get_var_env2('nomesh',undef);
$mkmesh = 0 if $nomesh;
$CASE = $case;
get_var_env2('N',4);
$nelem = $N**3.0;
$cs = 0.2*$nelem;		# desired chunk_size
$chunk_size = ($cs < 5000 ? 5000 : $cs > 40000 ? 40000 : $cs);
$hratio = 10;
$use_tetra = 1;
$maxits = 300;
$tol = 1e-7;
get_var_env2('subpart',350);
$subpart_entered = $subpart;
get_var_env2('NP',1);
$subpart = ceil($subpart/$NP);
$iisd_subpart = $subpart;
$nu = 1e-2;
get_var_env2('nlay',1);
get_var_env2('isp_maxits',12);

transcript("", qw(viscosity Dt));	# print variables on output
@vars= qw(CASE U L N Re viscosity 
    h Co Dt hratio leaky tol use_prismatic use_tetra);
octave_export_vars(">data.m.tmp",@vars);
if ($mkmesh) { system "octave -qH mkmesh.m > mkmesh.output.tmp"; }

1;
