#!/usr/bin/perl

sub pushf {
    my ($file,$N) = @_;
    return unless -f $file;
    my $F = "$file.".($N-1);
    unlink $F if -f $F;
    for ($j=$N-1; $j>=0; $j--) {
	my $to = "$file.$j";
	my $from = "$file.".($j-1);
	rename $from,$to if -f $from;
    }
    rename $file,"$file.0";
}

foreach  $N (100,200,400) {
foreach  $subpart (0,1,2,4,8,16) {
for ($nlay=0; $nlay<=5; $nlay++) {
    printf("start nlay $nlay\n");
    system "make nlay=$nlay > job.out.tmp";
    open OUT,"job.out.tmp";
    my ($nlay,$befo);
    my @solving = ();
    while (<OUT>) {
	$nlay = $1 if /use_interface_full_preco_nlay -> (\d*)/;
	$befo = $1 if /Before solving linear system.*\[\S*\s*(\S*)\]/;
	$iter = $1 if /iteration (\S*) /;
	$maxits = $1 if /^\s*maxits .> (\d*)\s*$/;       
	if (/After solving linear system.*\[\S*\s*(\S*)\]/) { 
	    $after_sol = $1; 
	    if ($iter>=$maxits) {
		print "Seems that didn't converge iter $iter, maxits $maxits\n";
	    }
	    push @solving,"t",($after_sol-$befo),"iter",$iter; 
	}
	$after_comp = $1 if /After computing profile.*\[\S*\s*(\S*)\]/;
    }
    $total = $after_sol - $after_comp;
    close OUT;
    pushf "job.out.tmp",3;
    my @mess = "nlay",$nlay,@solving,'total',$total;
    open STAT,">>isp.stat";
    my $mess = join(' ',@mess)."\n";
    print STAT $mess;
    print "end $mess";
    close STAT;
}}}

