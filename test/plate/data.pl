$Re = 100;
$Co = 2;
$re_start = 0;

if (1) {   ## large mesh
    $Rint=1;
    $Rext=2;
    $L = 10;
    $Ntheta = 80;
    $Nr=20;
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
