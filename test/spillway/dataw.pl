$h = 1;                         # Mean water height
$Lx = 5;			# Length of domain (wave-length)
$Re = 100;			# Reynolds number
$Froude = 0.5;			# Froude number

$yratio = 4;			# refinement towards top and bottom
$ref = 1;                       # Refinement parameter
$Ny = 10*$ref;                  # Number of elements along vert. direction
$Nx = 50*$ref;                  # Number of elements along horiz. direction

$gravity = 9.8;                 # Gravity acceleration (magnitude)
$u_av = $Froude*sqrt($gravity*$h);
				# Mean velocity
$viscosity = $u_av*$h/$Re;	# Fluid viscosity

$slope = 3*$u_av*$viscosity/($gravity*$h**2);                 
				# simulated bottom slope. A term g_x = slope*g is added
				# adjusted so as to give the desired Re and Froude number

$patm = 0;                      # Pressure at FS
$p_const_bc = !$initia;		# Impose = patm on the free surface
$fs_relax = 1;                  # relaxation factor for the 
                                # free surface evolution eq. 

1;
