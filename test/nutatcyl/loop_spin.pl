#!/usr/bin/perl 
#$Id: loop_spin.pl,v 1.1 2002/08/18 00:27:59 mstorti Exp $

# Traces the curve Mz(Omega_spin) at
# Omega_nut = 620rpm, nut_angle=15deg,

$n = 2;				#  number of points in curve

#================================================================
$allforce = "cylinder.all_force_spin.tmp";

@W = ();
require '../../tools/math.pl';

$force = "cylinder.force.tmp";
open ALLF, ">$allforce" || die "can't open $allforce\n";
$r=0;
while(1) {
    $Omega = refine2(1000,9000,\@W);
    for $nuta (0,15) {
	if (-f $force) { unlink $force; }
	system "make Omega=$Omega Omega_nut=620 nu=0.0572 nuta=$nuta mesh=10-20-I run\n";
	open FORCE, $force || die "couldn't find $force\n";
	my $line = <FORCE>;
	close FORCE;
	print ALLF "$Omega $nuta $line";
    }
    last if ++$r > $n;
}
close ALLF;
