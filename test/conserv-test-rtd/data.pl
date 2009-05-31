# $Id: cubcav.pl,v 1.8 2006/01/27 20:51:51 mstorti Exp $
require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";# Initializes ePerl 
require "math.pl";

$ndim = 2;

## Fine mesh
$ref = 1;
$n = 20*$ref;
$nz = 10*$ref;
$radius = 0.05;
$Lz = $radius;
$rratio = 1;

$conserv_test = 0;
$hz_fs = 1;
$hz_bot = 1;
$Omega = 0;
$Re = 100;
$Co = 2;
$weak_form = 0;

$Vmean = 1;
$Vmax = 2*$Vmean;
$nu = 2*$radius*$Vmean/$Re;
$dpdz = $Vmean**2.0*16.0/($radius*$Re);
$dpdz *= 3.0/8.0 if $ndim==2;
$h = $Lz/$nz;
$Dt = $Co*$h/$Vmax;

$T = 20*$Lz/$Vmean;
$nstep = $T/$Dt;
$f0 = 0;
$f1 = $Lz;
$nints = 10;                    # nbr of intervals (not values!!)
$nel = int(2**$ndim);
$nelb = $nel/2;

sub fvals {
    my ($f0,$f1,$nints) = @_;
    print "# fvals args: @_\n";
    my $fvals = "";
    for (my $j=0; $j<=$nints; $j++) {
        my $val = $f0 + ($f1-$f0)/$nints*$j;
        $fvals .= "$val ";
    }
    return $fvals;
}

## positions of planes must be avoided to
## coincide with an inter-element plane
my $tol = 1e-5;
$vals = fvals($tol,$Lz-2.3*$tol,$nints);

@vars= qw(n ndim nel nelb rratio  radius  nz  Lz  hz_fs  hz_bot
          Omega  Re Vmean nu dpdz h Dt T nstep conserv_test);
octave_export_vars(">data.m.tmp",@vars);
transcript2(@vars);	# print variables on output

if ($mkmesh) { 
    system "octave -qH mkcyl2d.m > mkmesh.output.tmp"; }

1;
