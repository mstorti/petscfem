$Re = 100;
$Co = 2;

if (1) {   ## large mesh
    $Rint=1;
    $Rext=4;
    $L = 8;
    $Ntheta = 80;
    $Nr=40;
    $Nx=80;
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
$Rext2=2*$Rext;
$Next = $Nr/2;
$hmin = 2*$PI*$Rint/$Ntheta;

1;
