#! /usr/bin/perl -w
use strict;
my $nmax=7;

sub indx {
    my ($i,$j,$n,$m) = @_;
    return $i*$m+$j;
}

open GEMMCODE,">gemmcode.cpp";
print GEMMCODE "#define NMAX $nmax\n";

my @loads;
for (my $n=1; $n<=$nmax; $n++) {
    for (my $m=1; $m<=$nmax; $m++) {
        for (my $p=1; $p<=$nmax; $p++) {
            print GEMMCODE "//","-" x 80,"\n";
            print GEMMCODE "static void p_${n}_${m}_${p}(double *a,double *b,double *c) {\n";
            for (my $j=0; $j<$n; $j++) {
                for (my $l=0; $l<$p; $l++) {
                    printf GEMMCODE "c[%d] = ",indx($j,$l,$n,$p);
                    for (my $k=0; $k<$m; $k++) {
                        printf GEMMCODE "a[%d]*b[%d]",
                        indx($j,$k,$n,$m),
                        indx($k,$l,$m,$p);
                        print GEMMCODE "+" if $k!=$m-1;
                    }
                    print GEMMCODE ";\n"
                }
            }
            print GEMMCODE "}\n";
            push @loads,"gemm_fun_table_load($n,$m,$p,&p_${n}_${m}_${p});\n"
        }
    }
}
close GEMMCODE;

open LOADFUNS,">loadfuns.cpp";
# print LOADFUNS "\n\n//","-" x 80,"\n";
print LOADFUNS @loads;
close LOADFUNS;
