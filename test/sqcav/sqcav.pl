require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
get_var_env('weak_form',0);
get_var_env('case',undef);
get_var_env('meas_perf',undef);
$use_triangles=0;
$save_file = ($case ? "sqcav.$case.tmp" : "sqcav.state.tmp");
if ($case =~ /^fractional_step_re(\d*)$/) {
    $Re = $1;
    $case = 'weak_form_1';
}
## Test the 'G_body' feature. 
$g_body= 0;
if ($case eq 'g_body') {
    $case='lu';
    $g_body=1;
}
$U = 1;				# typical velocity scale
$L = 1;				# Side of the cavity
$N = 20;
$hratio = 4;			# refinement (h_center/h_wall)
$Re = 400 unless defined $Re;			# Reynolds number
$viscosity = $L/$Re;		# viscosity
$Dt = 0.05;			# time step
$solver = "petsc";
$preco_type = "jacobi";
$maxits = 20;
$nstep = 500;
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
$fractional_step = 1;
# if ($fractional_step) { $solver = "iisd"; }
#
#__END_TRANSCRIPT__
@vars = qw(U L N Re viscosity h Co Dt hratio g_body use_triangles);
transcript("", @vars);	# print variables on output
octave_export_vars(">data.m.tmp",@vars);

1;
