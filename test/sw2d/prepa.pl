#!/usr/bin/perl

if (1){
    open IN,"dx/sw2d_tmpl.dx";
    $in = join "",(<IN>);
    close IN;
#
    $nstates=500;
    $ini=0;$stop=$nstates;
    @nn_ele = ();
    for (my $k=$ini; $k<=$stop; $k++) { 
	$out = $in;
	$out =~ s/%NUMERO%/$k/ ;
	open OUT,">dx/plano_$k.dx";
	print OUT $out;
	close OUT;
    }
}

