#__TRANSCRIPT__

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
get_var_env('weak_form',0);
get_var_env('case',undef);
print "# case entering sqcav.pl: $case\n";
get_var_env('meas_perf',undef);
get_var_env('use_triangles',0);
$use_triangles = 0;

$press_term=0;
if ($case eq 'press_term') {
    $press_term=1;
    $case='weak_form_0';
}

$u_rini = 0.5;
if ($case eq 'read_ini') {
    $read_ini_case = 1;
    $case = 'iisd_sbp2';
    $save_file = "sqcav.read_ini.state.tmp";
}

$dx = 0;
if ($case =~ /^dx_(.*)/) {
    $dx = 1;
    $dx_case = $1;
    $case = "weak_form_0"; 
    ${$dx_case} = 1;
    $synchro = 1 if $allf;
    $solver = 'iisd';
    $preco_type = 'jacobi';
    $steady = 0;
    $nstep = 20;
}
$save_file = ($case ? "sqcav.$case.tmp" : "sqcav.state.tmp") unless $save_file;

if ($case eq 'other') {
    $case = 'weak_form_1';
    $Re = 1000;
    $hratio = 4;
    $N = 40;
}

if ($case eq 'zwproc2_ref') { $case = 'zwproc2'; }

if ($case eq 'zwproc2') {
    $case = 'zwproc';
    $iisd_subpart = 2;
}

if ($case eq 'zwproc_ref') { $case = 'zwproc'; }

if ($case eq 'zwproc') {
    $case = 'weak_form_1';
    $preco_type = 'jacobi';
    $N=10;
    $Re=1;
    $nstep = 5;
    $solver = "iisd";
    $iisd_subpart = 1 unless $iisd_subpart;
    $steady = 0;
    $maxits = 100;
}

$steady = 1 unless defined $steady;

## Test the 'G_body' feature. 
$g_body= 0;
if ($case eq 'g_body') {
    $case='lu';
    $g_body=1;
}
$U = 1;				# typical velocity scale
$L = 1;				# Side of the cavity
$N = 20 unless $N;
$hratio = 10 unless $hratio;			# refinement (h_center/h_wall)
$Re = 1000 unless $Re;		# Reynolds number
$viscosity = $L/$Re;		# viscosity
$Dt = 0.01;			# time step
$solver = "petsc" unless $solver;
$preco_type = "lu" unless $preco_type;
$maxits = 1 unless $maxits;
$nstep = 100 unless $nstep;
#
$weak_form = 0;
$cache_gdu = 0;
$update_jacobian =0;
# `scalar' avoids returning the first match, which can be zero
# and then can be interpreted as `false' 
if (scalar $case =~ /weak_form_(\d)/) {
    $weak_form = $1;
} elsif ($case eq 'lu') {
    $solver = "petsc";
    $preco_type = "lu";
    $maxits = 1;
    $nstep = 2;
} elsif ($case =~ /^iisd/) {
    $solver = "iisd";
    $iisd_subpart = 2;
    $maxits = 200;
    $nstep = 2;
    $iisd_subpart = $1 if $case=~/_sbp(\d*)/;
    if ($case=~/_uj/) {
	$update_jacobian = 1;
	$nstep = 2000;
	$maxits=100;
	$iisd_subpart = 20;
	$N = 200;
    }
    $preco_type = "jacobi";
} elsif ($case eq 'big') {
    $N = 40;
#      $solver = "iisd";
#      $preco_type = "jacobi";
    $solver = "petsc";
    $preco_type = "lu";
    $maxits = 1;
    $nstep = 20;
    $cache_gdu = 0;
    $iisd_subpart = 10;
} else {
    die "don't know case $case\n";
}
$use_qharmm = 0;
$nstep = 0 if defined $read_ini_case;
#
#__END_TRANSCRIPT__
@vars = qw(U L N Re viscosity h Co Dt hratio g_body use_triangles u_rini);
transcript("", @vars);	# print variables on output
octave_export_vars(">data.m.tmp",@vars);

1;
