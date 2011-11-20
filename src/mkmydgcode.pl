#! /usr/bin/perl -w
use strict;
my $nmax=10;

sub indx {
    my ($i,$j,$n,$m) = @_;
    return $i*$m+$j;
}

my $warning = <<'EOM';
// This file was automatically created by "mkcode.pl"
// It is supposed not to be modified by hand. 
// It is not a standard header file. You should
// not include it in a program, this is only
// included once in file "fm2prod2.cpp"
EOM

open GEMMCODE,">mygmcode.h";
print GEMMCODE $warning;

open DEFFUNS,">mygmdefs.h";
print DEFFUNS $warning;
print DEFFUNS "#define PF_MYDGEMM_NMAX $nmax\n";
my @loads;
for (my $n=1; $n<=$nmax; $n++) {
    for (my $m=1; $m<=$nmax; $m++) {
        for (my $p=1; $p<=$nmax; $p++) {
            print GEMMCODE "//","-" x 80,"\n";
            my $fun = "p_${n}_${m}_${p}";
            ## print GEMMCODE "void prod2_subcache_t::$fun(double *a,double *b,double *c) {\n";
            print GEMMCODE "DEFFUN2($fun) {\n";
            print DEFFUNS "DECLFUN($fun);\n";
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
            ## push @loads,"gemm_fun_table_load($n,$m,$p,&p_${n}_${m}_${p});\n"
            push @loads,"LOADFUN($n,$m,$p);\n"
        }
    }
}
close GEMMCODE;
close DEFFUNS;

open LOADFUNS,">mygmload.h";
print LOADFUNS $warning;
# print LOADFUNS "\n\n//","-" x 80,"\n";
print LOADFUNS @loads;
close LOADFUNS;
