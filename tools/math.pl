# $Id: math.pl,v 1.6 2002/08/18 15:03:51 mstorti Exp $
# Some useful definitions for use with ePerl
#

$PI=2*atan2(1.,0.);
$E=exp(1.);

sub round {
    my $x = shift();
    my $y = shift();
    $y = 1 unless $y;
    return floor($x/$y + 0.5);
}

sub floor {
    my $x = shift();
    my $y = shift();
    $y = 1 unless $y;
    my $rr = $x/$y;
    my $r = int($rr);
    $r-- if $rr<0;
    return $r;
}

sub ceil {
    return floor(@_)+1;
}

sub tan { return sin($_[0])/cos($_[0]); }

sub asin {
    my $x=shift();
    return atan2($x,sqrt(1-$x*$x));
}
  
sub acos {
    my $x=shift();
    return atan2(sqrt(1-$x*$x),$x);
}

sub sinh {
    my $x=shift();
    return (exp($x)-exp(-$x))/2.;
}
  
sub cosh {
    my $x=shift();
    return (exp($x)+exp(-$x))/2.;
}

sub tanh {
    my $x=shift();
    return sinh($x)/cosh($x);
}
  
sub acosh {
    my $y=shift();
    die "acosh: argument must be >1!\n" if $y<1;
    return log($y+sqrt($y*$y-1));
}
  
sub asinh {
    my $y=shift();
    return log($y+sqrt($y*$y+1.));
}

sub atanh {
    my $y=shift();
    die "atanh: argument x must be |x| < 1.\n" if abs($y)>=1;
    return log((1+$y)/(1-$y))/2.;
}


sub asin {
    my $y=shift();
    return atan2($y,sqrt(1-$y**2));
}

sub log10 {
    my $y=shift();
    return log($y)/log(10.);
}

# Returns a number floored to a grid of the form 1,2,5,10,20,50, etc...
# The grid may be changed.
#
# Usage: log_floor($x)           -   floors to the standard grid 1,2,5,10,...
#        log_floor($x,'fine')    -   floors to a fine grid 1,1.5,2,3,4,5,7,10,...
#        log_floor($x,[1,1.5,2,3,4,5,6,7,8,9,10])  -   floors to a user defined  grid 
#                                                      the grid starts in 1 and ends in 10. 
sub log_floor {
    my $x = shift();
    return $x if $x==0;
    my $sign=1;
    if ($x<0) {
	$sign=-1;
	$x = -$x;
    }
    $log10 = log($x)/log(10.);
    $expo = floor($log10);
    $man = 10**($log10-$expo);
    my $grid = shift();
    $grid = [1,1.5,2,3,4,5,7,10] if $grid eq 'fine';
    $grid = [1,2,4,8,10] if $grid eq 'binary';
    $grid = [1,2,5,10] unless $grid; # default grid
    for($j=1; $j<=$#{$grid}; $j++) {
	$r = $grid->[$j];
	return $sign*$grid->[$j-1]*(10**$expo) if $man<$r;
    }
}

# Returns a number ceiled to a grid of the form 1,2,5,10,20,50, etc...
# The grid may be changed.
#
# Usage: see log_floor
#
sub log_ceil {
    my $x = shift();
    return $x if $x==0;
    my $sign=1;
    if ($x<0) {
	$sign=-1;
	$x = -$x;
    }
    $log10 = log($x)/log(10.);
    $expo = floor($log10);
    $man = 10**($log10-$expo);
    my $grid = shift();
    $grid = [1,1.5,2,3,4,5,7,10] if $grid eq 'fine';
    $grid = [1,2,4,8,10] if $grid eq 'binary';
    $grid = [1,2,5,10] unless $grid; # default grid
    foreach $r (@{$grid}) {
	return $sign*$r*(10**$expo) if $man<$r;
    }
}

sub refine {
    $nup=shift();
    my @nu=@{$nup};
    my $nu1=shift();
    my $nu2=shift();
    if ($#nu==-1) {
	$nu=$nu1;
#	print "paso por 1\n";
    } elsif ($#nu==0) {
	$nu=$nu2;
#	print "paso por 2\n";
    } else {
	for ($k=0;$k<$#nu;$k++) {
#	    print "k: $k\n";
	    $dif=$nu[$k+1]-$nu[$k];
#	    print "dif: $dif\n";
	    if ($k==0 || $dif>$difmax) {
		$difmax=$dif;
		$kmin=$k;
	    }
	}
#	print "intervalo minimo entre ",$nu[$kmin]," y ",$nu[$kmin+1],"\n";
	$nu=($nu[$kmin]+$nu[$kmin+1])/2;
    }
    push @nu,$nu;
    @nu = sort @nu;
    @{$nup}=@nu;
    return $nu;
}

sub max {
    my $maxx = $_[0];
    for (my $j=1; $j<=$#_; $j++ ) {
	$maxx = $_[$j] if $_[$j]>$maxx;
    }
    return $maxx;
}

sub min {
    my $minn = $_[0];
    for (my $j=1; $j<=$#_; $j++ ) {
	$minn = $_[$j] if $_[$j]<$minn;
    }
    return $minn;
}

sub irand {
    if ($#_ == 0) { return int(rand()*$_[0]); }
    elsif ($#_ == 1) { return $_[0]+irand($_[1]-$_[0]+1); }
    else { die "usage: irand(n) or irand(i,j)\n"; }
}

# usage: ($dist,$nearer) = distance($x,\@v)
# compute nearer point from values in @v to $x
sub distance {
    my ($x,$v) = @_;
    my $dist_min, $nearer, $dist;
    for ($k=0; $k<=$#$v; $k++) {
	$dist = abs($x-$v->[$k]);
	if (!defined $dist_min || $dist < $dist_min) {
	    $dist_min = $dist;
	    $nearer = $v->[$k];
	}
    }
    return $dist_min,$nearer;
}

=head1 NAME

refine2: finds a new point where to compute a value

=head1 SYNOPSIS

Usage: C<@v_ref = refine2($v1,$v2,\@v); >

Finds the next point to compute a value in the range C<$v1,$v2>, if
you already computed at @v. It finds the point in the interval
C<[$v1,$v2]> that is at the larger distance from all points in C<@v>. 

=head1 DESCRIPTION

Suppose you are running a program for a series of values of parameter
C<x> and recover the value of magitude C<y>. You already computed the
values in C<@v> and want to compute a new value. C<refine2> finds the
best point in the interval C<[$v1,$v2]> where to compute the next
value, that is the point int the interval that is at the larger
distance from all points in C<@v>. The algorithm is based in noting
that the point should be either one of the extremes of the interval or
either the midpoints of those points in C<@v> that fall in the interval. 
For instance if you call 

    @v = ();
    while (...) { @v = refine2(0,1,\@v; }

then the sequence generated is C<{0,1,0.5,0.25,0.75,0.125,0.375,...}>.
Restart: Note that if the points in C<@v> may fall outside the
interval. This is useful if you computed at a some values, then have
some interest at some interval, for instance C<0.3,0.5]>, then you
simply restart with 

    while (...) { @v = refine2(0.3,0.5,\@v; }

and this will find optimal points in the given interval. 

=head1 AUTHOR

Mario A. Storti C<mstorti@intec.unl.edu.ar>

=cut

sub refine2 {
    my ($v1,$v2,$v) = @_;
    my $tol = 1e-6*($v2-$v1);
    my @v = @$v;
    my $newv;
    # eliminate all values outside [$v1,$v2]
    my @vv = ();
    foreach $w (@v) { push @vv,$w if $w > $v1-$tol && $w < $v2+$tol; }

    # put tentative evaluation points
    @v = ();
    push @v,$v1,$v2;
    for (my $j=0; $j<$#vv; $j++) { push @v,($vv[$j]+$vv[$j+1])/2.; }
    undef @vv;
    my $dmax, $wref;
    foreach $w (@v) {
	my ($d,$n) = distance($w,$v);
	if (!defined $dmax || $d > $dmax) {
	    $wref = $w;
	    $dmax = $d;
	}
    }
    my @w = @$v;
    push @w,$wref;
    @$v = sort @w;
    return $wref;
}

# other version of refine.
sub refine3 {
    my ($v1,$v2,$v) = @_;
    my @v = @$v;
    my $newv;
    if ($#v==-1) {
	$newv = $v1;
    } elsif ($#v==0) {
	$newv = $v2;
    } else {
	my $maxdif = 0;
	for (my $k=0; $k<$#v; $k++) {
	    my $dv = $v[$k+1]-$v[$k];
	    if ($dv > $maxdif) {
		$maxdif = $dv;
		$newv = ($v[$k+1]+$v[$k])/2;
	    }
	}
    }
    push @$v,$newv;
    @$v = sort {$a <=> $b } @$v;
    return $newv;
}

1;
