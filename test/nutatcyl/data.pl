#$Id: data.pl,v 1.1 2002/08/17 22:27:47 mstorti Exp $
$omega = sqrt($PI);             # Typical sloshing frequency
$Omega = 2.;			# Rotation velocity
$nstep = 2000;			# Number of time steps
$Dt = 0.1;			# Time step
$tend = $nstep * $Dt;		# Computation stops here
$TT=0.2;			# Impulse ends here
$verbose = 0;
1;
