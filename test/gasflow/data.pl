$flat_plate = 0;
$slip = 1;

## GEOMETRY
$L = 8;				# semi-length of nozzle
$Rin = 1;			# radius at nozzle inlet/outle
$Rn = 0.5;			# radius at throttle (minimum radius)
$nw = 3;			# throttle half-width. The shape is a `cosh' and
				# at this distance from the center it recovers 90% of the
				# inlet width.
$Nr = 5;			# number of elements in radius direction
$Nx = 40;			# number of elements in x direction
$r_ratio = 8;			# hr(axis)/hr(wall)
$x_ratio = 16;			# hx(inlet) / hx(throttle)

# BOUNDARY CONDITIONS
$Mach_ref   = 0.4;
$rho_ini    = 1.15;
$p_ini      = 1.0e5; 
$p_ref      = 1.0e5;
$rho_ref    = 1.15;

$cc_ref     = sqrt($ga*$p_ref/$rho_ref);
$u_ref      = $Mach_ref*$cc_ref;
$u_ini      = $u_ref;
$v_ini      = 0;
$v_ref      = $v_ini;

$compute_s = 1;
$compute_h = 1;

1;
