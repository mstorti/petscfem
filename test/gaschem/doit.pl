#!/usr/bin/perl

@procs = qw(node17 node14 node15 node18 node19 
	    node20 node21 node22 node23);

for ($j=0; $j<@procs; $j++) {
    $noproc = $procs[$j];
    print "---<*>" x 10,"\n";
    print "---<*>" x 10,"\n";
    print "---<*>" x 10,"\n";
    print "doit.pl: lanzando sin proc $noproc\n";
    open PROCT,">proctable";
    print PROCT "node1 1.0 server\n";
    for ($k=0; $k<@procs; $k++) {
	if ($k!=$j) { print PROCT "$procs[$k] 1.0\n"; }
    }
    system "/usr/bin/make run | tee nohup.log";
#    system "echo \"---- proctable: -----\" >> nohup.log"; 
    system "cat proctable >> nohup.log";
    rename "nohup.log","nohup.log.noproc-$noproc";
}
