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

# foreach  $N (100,200,400) {
# for ($nlay=0; $nlay<=5; $nlay++) {
do {{ $N=200; $nlay=1;
       foreach  $subpart (1,2,4,8,16) {
       my @mess = ("nlay",$nlay,"subpart",$subpart,"N",$N);
       my $mess = join(' ',@mess)."\n";
       printf("start $mess");
       system "make nlay=$nlay N=$N subpart=$subpart > job.out.tmp";
       open OUT,"job.out.tmp";
       my ($nlay,$befo);
       my @solving = ();
       my $iter = 0;
       my $maxmem = 0;
       while (<OUT>) {
	   $nlay = $1 if /use_interface_full_preco_nlay -> (\d*)/;
	   $befo = $1 if /Before solving linear system.*\[\S*\s*(\S*)\]/;
	   $iter = $1 if /iteration (\S*) /;
	   $maxits = $1 if /^\s*maxits .> (\d*)\s*$/;       
	   $maxmem = $1 if /After solving linear system.*\[Memory usage.*max (\S*),/;
	   $after_sol = $1 if (/After solving linear system.*\[\S*\s*(\S*)\]/);
	   if (/After solving linear system.*\[Memory usage.*max (\S*),/) { 
	       if ($iter>=$maxits) {
		   print "Seems that didn't converge iter $iter, maxits $maxits\n";
	       }
	       push @solving,"t",($after_sol-$befo),"iter",
	       $iter,"maxmem",$maxmem; 
	   }
	   $after_comp = $1 if /After computing profile.*\[\S*\s*(\S*)\]/;
       }
       $total = $after_sol - $after_comp;
       close OUT;
       pushf "job.out.tmp",3;

       push @mess,@solving,'total',$total; 
       $mess = join(' ',@mess)."\n";
       open STAT,">>isp.stat";
       print STAT $mess;
       print "end $mess";
       close STAT;
}}}

