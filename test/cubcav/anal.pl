#!/usr/bin/perl -w
#$Id: anal.pl,v 1.1 2002/09/16 19:34:31 mstorti Exp $
use strict;

my @time = ();
my $tprev;
%time=();
%steps=();
while (<>) {
    if (/^--.*After(.*)\s*--(.*)\s(\S*)\]/) {
	my $job1 = $1;
	my $t1 = $2;
	my $line1 = $_;
	while (<>) {
	    if ($tprev) { 
		if (/^--.*After(.*)\s*--(.*)\s(\S*)\]/) {
		    my $job2 = $1;
		    my $t2 = $2;
		    my $line2 = $_;
		    die "\"After\" line doesn't match \"Before\" line\n",
		    $line1,$line2 unless $job2 eq $job1;

		    $time{$job1} = 0 unless $time{$job1};
		    $time{$job1} += $t2-$t1;

		    $steps{$job1} = 0 unless $steps{$job1};
		    $steps{$job1}++;
		}
	    }
	}
    }

    while (my ($job,$time) = each %time) {
    }
