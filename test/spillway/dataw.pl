$h = 1;                         # Mean water height
$Lx = 5;			# Length of domain (wave-length)
$Re = 500;			# Reynolds number
$Froude = 2;			# Froude number

$yratio = 4;			# refinement towards top and bottom
$ref = 2;                       # Refinement parameter
$Ny = 10*$ref;                  # Number of elements along vert. direction
$Nx = 50*$ref;                  # Number of elements along horiz. direction

$gravity = 9.8;                 # Gravity acceleration (magnitude)
$u_av = $Froude*sqrt($gravity*$h);
$u_max = 1.5*$u_av;
				# Mean velocity
$viscosity = $u_av*$h/$Re;	# Fluid viscosity

$slope = 3*$u_av*$viscosity/($gravity*$h**2);                 
				# simulated bottom slope. A term g_x = slope*g is added
				# adjusted so as to give the desired Re and Froude number

$patm = 0;                      # Pressure at FS
$p_const_bc = !$initia;		# Impose = patm on the free surface
$fs_relax = 1;                  # relaxation factor for the 
                                # free surface evolution eq. 

$Dt = 0.05;			# time step
$eta0 = 0.05;			# amplitude of free surface elevation
				# perturbation
$dx_steps = 1;			# visualize with DX

$L_bump = 1;			# length of the bump
$t = 0.1;			# thickness
$restart = 0;			# flags whether to restart a previous run

1;
