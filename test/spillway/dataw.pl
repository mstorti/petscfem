$h = 1;                         # Mean water height
$Lx = 20;			# Length of domain (wave-length)
$Re = 1e2;			# Reynolds number
$Froude = 1.0;			# Froude number

$yratio = 4;			# refinement towards top and bottom
$ref = 1;                       # Refinement parameter
$Ny = 14*$ref;                  # Number of elements along vert. direction
$Nx = 200*$ref;			# Number of elements along horiz. direction

$gravity = 9.8;                 # Gravity acceleration (magnitude)
$u_av = $Froude*sqrt($gravity*$h);
$u_max = 1.5*$u_av;
				# Mean velocity
$viscosity = $u_av*$h/$Re;	# Fluid viscosity

$slope = (3*$u_av*$viscosity/($gravity*$h**2));                 
				# simulated bottom slope. A term g_x = slope*g is added
				# adjusted so as to give the desired Re and Froude number

$patm = 0;                      # Pressure at FS
$p_const_bc = !$initia;		# Impose = patm on the free surface
$fs_relax = 1.0;		# relaxation factor for the 
                                # free surface evolution eq. 
$fs_smoothing_coef = 0.1;	# FS smoothing coeficient (heat equation like)
$nfilt = 20;			# Number of times the filter is applied

$Dt = 0.01;			# time step
$eta0 = 0.05;			# amplitude of free surface elevation
				# perturbation
$dx_steps = 0;			# visualize with DX

$L_bump = 2;			# length of the bump
$t = 0.2;			# thickness

$restart_dir = "RUN";
$restart_step = -1;		# restart from a previous run

$restart = $restart_step>=0;	# flags whether to 

1;
