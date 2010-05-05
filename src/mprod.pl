#! /usr/bin/perl -w
use strict;
use File::Slurp;

my $prodtmpl = "FastMat2 & prod(\n__MATS__const int m,INT_VAR_ARGS);\n";
my $mprodtmpl = read_file('./mprodtmpl.h');

open MPHDR,">./mprod.h";
open MPHDR2,">./mproddef.h";
my $nmatmax=10;
for (my $nmat=3;$nmat<=$nmatmax; $nmat++) {
    my $mats = "";
    my $pack = "";
    for (my $k=0; $k<$nmat; $k++) {
        $mats .= "                const FastMat2 & A$k,\n";
        $pack .= "  mat_list.push_back(&A$k);\n"
    }
    my $t = $prodtmpl;
    $t =~ s/__MATS__/$mats/;
    print MPHDR "\n\n",$t;

    $t = $mprodtmpl;
    $t =~ s/__MATS__/$mats/;
    $t =~ s/__PACK__/$pack/;
    print MPHDR2 "\n\n",$t;

}
close MPHDR2;
close MPHDR;
