#
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
    for (my $j=1; $j<$#_; $j++ ) {
	$maxx = $_[$j] if $_[$j]>$maxx;
    }
    return $maxx;
}

sub min {
    my $minn = $_[0];
    for (my $j=1; $j<$#_; $j++ ) {
	print "minn: cur val: $_[$j]\n";
	$minn = $_[$j] if $_[$j]<$minn;
    }
    return $minn;
}

1;
