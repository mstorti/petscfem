#!/usr/bin/perl

sub averg {
    my $sum=0;
    my $count = 0;
    my $ref = shift();
    foreach my $t (@$ref) {
	$sum += $t;
	$count++;
    }
    return $sum/$count;
}

#open IN,"isp30c.stat";
#$stream = IN;

$stream = STDIN;
print "nlay sbp ispits mem iter av. tav\n";
while (<$stream>) {
    my @items = split(" ",$_);
    my %items;
    while (@items) {
	my $key = shift @items;
	my $val = shift @items;
	if (!exists($items{$key})) { $items{$key} = $val; }
	else {
	    if (!ref($items{$key})) { $items{$key} = [$items{$key}]; }
	    push @{$items{$key}}, $val;
	}
    }
    printf "%2d %3d %3d %7.1f %7.2f %7.1f\n",$items{nlay},$items{subpart},
    $items{isp_maxits},averg($items{maxmem})/1000.0,averg($items{iter}),$items{tav};
}
