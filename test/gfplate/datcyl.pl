#! /usr/bin/perl

require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$R = 1;				# Radius of cylinder
$Rext = 5;			# radius of external cylinder
				# (external boundary)

$Nr = 20;			# Number of nodes along radial direction
$Nphi = 2*$Nr;			# Number of nodes along skin
$rratio = 5;			# Refinement in radial direction

$Machin = 0.5;			# Mach at inlet
$gamma = 1.4;			# Cp/Cv for gas
$Rgas = 1;			# Gas constant
$rhoref = 1;			# Reference density
$Tref = 1;			# Reference temperature
$Twall = 1.5;			# Temperature at wall

$abso = 1;			# Use absorbing b.c.'s?
$restart = 0;			# Is this a restart?
$use_symm = 0;			# Use symmery?
$use_twall = 0;			# Impose T on cylinder skin?
$use_non_slip = 0;		# Impose v=0 on skin?
$dv_pert_symm = 0.; 		# Non symmetric perturbation
				# to initiate unsteady flow

$pref = $rhoref*$Rgas*$Tref;
$cref = sqrt($gamma*$pref/$rhoref);
$uref = $Machin*$cref;

@vars = qw(R Rext Nr Nphi rratio Machin 
	   gamma Rgas rhoref Tref abso
	   pref cref uref use_symm use_twall
	   use_non_slip dv_pert_symm Twall);
octave_export_vars(">data.m.tmp",@vars);

$oscript = ($use_symm ? 'mkcyl' : 'mkcyl2');
system "octave -qH $oscript.m > $oscript.log.tmp";

if (!$dx && !$restart) { 
    system "echo -n > cylabso.some-rslt.tmp"; 
}

1;
