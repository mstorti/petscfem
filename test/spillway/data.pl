$C = 0.053547;			# Constants of spillway expr.
$E = 1.85;
$H1 = 66.50;			# vertical position of spillway top position
$L1 = 18;			# Distance from top of spillway to start
				# of flat bottom
$L2 = 50;			# Flat bottom length
$h1 = 4;                        # Water height at top of spillway
$y2 = 65.79;			# restitution height

$ref = 2;                       # refinement parameter
$Ny = 10*$ref;                  # Number of elements along vert. direction
$Nx = 50*$ref;                  # Number of elements along horiz. direction

$gravity = 9.8;                 # gravity acceleration (magnitude)
$uin = sqrt($gravity*$h1);      # velocity at inlet
$Dt = 0.1;
$patm = 0.;			# atmosferix pressure (taken at the
                                # restitution height). Should take into
                                # account the velocity of the fluid at the outlet. 
$fs_relax = 1.;			# relaxation factor for the 
                                # free surface evolution eq. 
$x_damp = $L1+0.5*$L2;          # a wave damper term is added for `x > x_damp'
$initia = 0;			# run steady initialization job
#$steady_state_file = 
#    "./RUN1/spillway.state_21.tmp";# Read pressure for FS nodes from this file
$p_const_bc = !$initia;		# Impose = patm on the free surface
$fs_relax = 0 if !$p_const_bc;

1;
