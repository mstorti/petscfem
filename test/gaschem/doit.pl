#!/usr/bin/perl

@procs = qw(node17 node14 node15 node18 node19 
	    node20 node21 node22 node23);

$j = shift();
print "doit.pl: lanzando sin proc $procs[$j]\n";
open PROCT,">proctable";
print PROCT "1.0 node1 server\n";
for ($k=0; $k<@procs; $k++) {
    if ($k!=$j) { print PROCT "1.0 $procs[$k]\n"; }
}
print "doit.pl: done\n";
