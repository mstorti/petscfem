# $Id: cubcav.pl,v 1.2 2003/11/25 21:49:10 mstorti Exp $
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
$subpart = ceil($subpart/$NP);
$iisd_subpart = $subpart;
$nu = 1e-2;
get_var_env2('nlay',1);
get_var_env2('isp_maxits',12);
$gaschem = 1;
if ($gaschem) {
    $gravity = 9.8;		# gravity acc. [m/sec2]
    $xO_in = 0.2;		# molar fraction O2 at inlet [1]
    $KO = 1.3516e-5;		# Henry constant O2 - N2 [1]
    $KN = 0.6788e-5; 		# 
    $Rgas = 8.314;		# Gas const    
    $Tgas = 293.0;		# Temperature    [K]
    $nu_t = 0.001;		# turbulent kinematic
				# viscosity   [m2/sec]
    $hm_fac = 1.0;		# factor controlling
				# gas-liquid exchange [1]
    $rb = 0.01;			# radius of bubble at inlet
    $alpha_in = 0.01;		# gas fraction
    $patm = 1e5;		# Atmosferic pressure
    $Dt = 0.1;			# 
    $sat_in = 0.3;		# Saturation at inlet

    $vb = 4/3*$PI*$rb**3.0;	# bubble volume
    $Nb_in = $alpha_in/$vb;	# Nbr of bubbles at inlet
    $C_in = $patm*$alpha_in/($Rgas*$Tgas); # 
				# Total conc. of gas at inlet

    $CO_in = $C_in*$xO_in;	# Conc. of O2 at inlet
    $CN_in = $C_in-$CO;		# Conc. of N2 at inlet
    $CdO_in = $sat_in*$KO*$patm*$xO_in;	# 
				# unperturbed conc. of diss O2 at inlet
    $CdN_in = $sat_in*$KN*$patm*(1.0-$xO_in); # 
				# unperturbed conc. of diss N2 at inlet

    $Nb_scale = 1e3*$Nb_in;	# Scales Nb for reducing difference
				# in magnitude between variables 
    $uref = 5;			# Velocity scale 
    $L = 5;			# Length scale
    $Nb = 2e4;			# Typical number of bubbles
    $Nb_source = $Nb*$uref/$L;	# Needed source

    $Nb_ctff = 0.1*$Nb_in;
    $CN_ctff = 0.1*$CN_in;
    $CO_ctff = 0.1*$CO_in;
    $CdO_ctff = 0.001*$CdO_in;
    $CdN_ctff = 0.001*$CdN_in;
    $N = 4;
}

1;
