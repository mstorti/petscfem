$Re = 1000;			# Reynolds number
$Co = 2;			# Courant number
$re_start = 0;			# restart previous run? 

$ev = sub { print $_[0]; };

$omega=0.25;			# filter param (time scale)
$xif=1;				# filter param (damping)
/`/;
$filter1 = <<EOM;
# Filter
n_coef 3
a_coef &$ev(1./$omega**2+$xif/$omega+0.5) \
       &$ev(-2/$omega**2)                 \
       &$ev(1./$omega**2-$xif/$omega+0.5)
b_coef &$ev(0.5+$xif/$omega) \
       0.                   \
       &$ev(0.5-$xif/$omega)
EOM

$omega=0.;			# filter param (time scale)
$filter2 = <<EOM;
# Filter
n_coef 2
a_coef 1 &$ev(1-$omega)
b_coef $omega 0
EOM
/`/;

$filter = $filter2;

if (1) {   ## large mesh
    $Rint=1;
    $Rext=3;
    $L = 10;
    $Ntheta = 80;
    $Nr=30;
    $Nx=100;
    $du_ini_pert = 0.3;
} else {  ## small mesh, small domain
    $Rint=1;
    $Rext=2;
    $L = 6;
    $Ntheta = 40;
    $Nr=20;
    $Nx=40;
    $du_ini_pert = 0.3;
}
$Rext2=10;
$Next = 3*$Nr/2;
$hmin = 2*$PI*$Rint/$Ntheta;
$Dt=$Co*$hmin;

1;
