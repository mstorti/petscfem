#$Id: vaughn.pl,v 1.1 2002/08/17 22:27:47 mstorti Exp $
$data = $ENV{'data'};
$data =~ s/\#/\$/g;
print "# data: $data\n";
eval($data);

get_var_from_env('Omega',9000);	# Spin velocity in r.p.m.
get_var_from_env('Omega_nut',
		 620);		# nutation velocity r.p.m.
get_var_from_env('nuta',15);	# nutation angle [degree]
get_var_from_env('nu',0.1);	# [m^2/sec] kinematic viscosity

# Convert to 1/sec
$Omega = $Omega/60.*2*$PI;	# Spin velocity [1/sec]
$Omega_nut 
    = $Omega_nut/60.*2*$PI;	# Nutation velocity [1/sec]

$nutation_angle 
    = $nuta*$PI/180.;		# Nutation angle
$nstep = 1;			# Number of time steps
$Dt = 0.2;			# Time step
$verbose = 0;

$RRin = 2.375;
$Lzin = 20.75;
$RR = $RRin*2.54/100.;		# [m] radius of cylinder
$cLz = $Lzin/$RRin;		# height of cylinder
$rho = 1;			# [kg/m^3] Density
