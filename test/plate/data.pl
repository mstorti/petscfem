$Re = 1000;			# Reynolds number
$Co = 2;			# Courant number
$re_start = 0;			# restart previous run? 
$dx_steps = 1;
$dx_read_state_from_file = 0;
$nsaverot = 1;

sub print_filter {
    my ($a,$b) = @_;
    my @r;
    push @r,"n_coef ",scalar @$a,"\na_coef ";
    my $aa,$bb;
    for $aa (@$a) { push @r, " $aa"; }
    push @r,"\nb_coef";
    for $bb (@$b) { push @r, " $bb"; }
    push @r,"\n";
    return join('',@r);
}

$omega=0.25;			# filter param (time scale)
$xif=1;				# filter param (damping)
@a = (1./$omega**2+$xif/$omega+0.5,-2/$omega**2,
      1./$omega**2-$xif/$omega+0.5);
@b = (0.5+$xif/$omega,0.,0.5-$xif/$omega);
$filter1 = print_filter(\@a,\@b);

$omega=0.;			# filter param (time scale)
@a = (1,-1+$omega);
@b = ($omega,0);
$filter2 = print_filter(\@a,\@b);

$filter = $filter1;

if (0) {   ## large mesh
    $Rint=1;
    $Rext=2;
    $L = 10;
    $Ntheta = 80;
    $Nr=30;
    $Nx=100;
} elsif (0) {   ## intermediate
    $Rint=1;
    $Rext=2;
    $L = 8;
    $Ntheta = 48;
    $Nr=20;
    $Nx=60;
} else {  ## small mesh, small domain
    $Rint=1;
    $Rext=2;
    $L = 4;
    $Ntheta = 24;
    $Nr=10;
    $Nx=16;
}
$du_ini_pert = 0.3;
$Rext2=10;
$Next = 3*$Nr/2;
$hmin = 2*$PI*$Rint/$Ntheta;
$Dt=$Co*$hmin;

1;
