#!/usr/bin/perl -w
use strict;

my ($step,$fifoin,$fifoout) = @ARGV;

open FIFOIN,">$fifoin";
print FIFOIN "$step\n";
close FIFOIN;

open FIFOOUT,"$fifoout";
my $ans = <FIFOOUT>;
# print "ans: $ans\n";
die "can't find OK on fifo!!\n" unless $ans =~ /^OK/;
close FIFOOUT;
