#!/usr/bin/perl

print "nlay sbp N  ispits mem  tav\n";
while (<STDIN>) {
    my @items = split(" ",$_);
    my %items;
    while (@items) {
	my $key = shift @items;
	my $val = shift @items;
	if (!exists($items{$key})) { $items{$key} = $val; }
	else {
	    if (!ref($items{$key})) { $items{$key} = [$items{$key}]; }
	    push @$items{$key}, $val;
	}
    }
    if (/^nlay (\d)* subpart (\d*) N (\d*) isp_maxits (\d*) iter (\d*) .* iter (\d*) .* maxmem (\d*) tav (\S*) total/) {
	$itav = 
	printf "%2d %3d %3d %3d %7.1f %7.1f\n",$1,$2,$3,$4,$5/1000.0,$6;
    }
}
