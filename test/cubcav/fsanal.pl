#!/usr/bin/perl -w
#$Id: fsanal.pl,v 1.2 2002/09/13 21:04:34 mstorti Exp $
use strict;

my @time = ();
my $tprev;
while (<>) {
    if (/^--.*After(.*)\s*--.*\s(\S*)\]/) {
	if ($tprev) { 
	    my $dt = $2-$tprev;
	    push @time,$dt; 
	    chomp;
#	    printf "%s -> %f\n",$_,$dt;
	}
	$tprev = $2;
    }  
}

for (my $j=0; $j<8; $j++) { shift @time; }

my @t = (0.0) x 7;
my @t2 = (0.0) x 7;
my $nt = 0;
while (1) {
    last if ($#time+1)<7;
    my $sum=0;
    for (my $j=0; $j<7; $j++) { 
	my $v = shift(@time); 
	printf "%f ",$v;
	$t[$j] += $v; 
	$t2[$j] += $v*$v; 
	$sum += $v;
    }
    $nt++;
    print " tot: $sum\n";
}
print "\nStat: -----------\n";
my $sum=0;
for (my $j=0; $j<7; $j++) { 
    $t[$j] /= $nt;
    $sum += $t[$j];
    $t2[$j] /= $nt;
    printf "%7.2f (sigma %7.2f)\n",$t[$j],($t2[$j]-$t[$j]**2.)/$nt; 
}
printf "total %7.2f (number of steps: $nt)\n",$sum;


