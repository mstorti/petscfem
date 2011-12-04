#! /usr/bin/perl -w
use strict;
use File::Slurp;

my $prodtmpl = "FastMat2 & prod(\n__MATS__const int m,INT_VAR_ARGS);\n";
my $mprodtmpl = read_file('./mprodtmpl.h');

my $warning = <<'EOM';
// This file was automatically created by mprod.pl
// It is supposed not to be modified by hand. 
// It is not a standard header file. You should
// not include it in a program, this is only
// included once in file fastmat2c.cpp
EOM

open MPHDR,">./mprod.h";
open MPHDR2,">./mproddef.h";
print MPHDR $warning;
print MPHDR2 $warning;

my $nmatmax=10;
for (my $nmat=2;$nmat<=$nmatmax; $nmat++) {
    my $mats = "";
    my $pack = "";
    my $ctx_check = "";
    for (my $k=0; $k<$nmat; $k++) {
        $mats .= "                const FastMat2 & A$k,\n";
        $pack .= "  mat_list.push_back(&A$k);\n";
        $ctx_check .= "    ctx->check(&A$k);\n";
    }
    my $t = $prodtmpl;
    $t =~ s/__MATS__/$mats/;
    print MPHDR "\n\n",$t;

    $t = $mprodtmpl;
    $t =~ s/__MATS__/$mats/;
    $t =~ s/__PACK__/$pack/;
    $t =~ s/__CTX_CHECK__/$ctx_check/;
    $t =~ s/__NMAT__/$nmat/;
    print MPHDR2 "\n\n",$t;

}
close MPHDR2;
close MPHDR;
