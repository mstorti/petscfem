# $Id: cubcav.pl,v 1.14.2.1 2004/05/22 02:48:27 mstorti Exp $
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
require "math.pl";

get_var_env2('case','laplace');
$case_in = $case;
$srfgath = 0;
if ($case eq 'srfgath') {
    $srfgath = 1;
    $case = 'plain';
}
$bupl = 1;
if ($case =~ /plain_bupl(\d)/) {
    $bupl = $1;
    $case = 'plain';
}
get_var_env2('data_dir','.');
$endur_test = $data_dir ne '.';
get_var_env2('job',undef);
get_var_env2('mkmesh',1);
get_var_env2('nomesh',undef);
$mkmesh = 0 if $nomesh;
$CASE = $case;
get_var_env2('N',30);
$nelem = $N**3.0;
$N = 14 if $case eq 'plain';
$cs = 0.2*$nelem;		# desired chunk_size
$chunk_size = ($cs < 5000 ? 5000 : $cs > 40000 ? 40000 : $cs);
$hratio = 3;
$geom = 'tetra';
$geom = 'cartesian3d' if $case eq 'plain';
$maxits = 300;
$tol = 1e-7;
get_var_env2('subpart',350);
$subpart_entered = $subpart;
get_var_env2('NP',1);
$subpart = ceil($subpart/$NP);
$iisd_subpart = $subpart;
get_var_env2('nlay',1);
get_var_env2('isp_maxits',12);
get_var_env2('nu',1/100.);
$nnwt = 5;
$elemset = 'nsi_tet_les_fm2';

if ($srfgath) {
    $N = 4;
    $hratio = 3;
    $geom = 'tetra';
    $nnwt = 1;
    $elemset = 'qharmm';
}

$use_tetra = $use_prismatic = 0;
$use_prismatic = 1 if $geom eq 'prismatic';
$use_tetra = 1 if $geom eq 'tetra';
$npg = ($use_tetra ? 4 : $use_prismatic ? 6 : 8);
$nel = $npg;

$icone = $data_dir.'/cubcav.con'.
    ($use_tetra ? '-tetra' 
     : $use_prismatic ? '-prism' : '').'.tmp';

@vars= qw(CASE U L N Re viscosity 
    h Co Dt hratio leaky tol use_prismatic use_tetra 
	  case_in case_oct npg nel icone elemset);
octave_export_vars(">data.m.tmp",@vars);
if ($mkmesh) { system "octave -qH mkmesh.m > mkmesh.output.tmp"; }


1;
