#!/usr/bin/perl
#$Id: split_mat.pl,v 1.2 2004/07/16 17:29:15 rodrigop Exp $
#$matvec = "lap_matrix10001000.m";
$matvec = "lap_matrix500500.m";
$mat = "lap_A.dat";
$vec = "lap_b.dat";
#lap 1000 x 1000
#$n = 8986009; #nonzeros
#$N = 1002001; #mat dim
$n = 2243009; #nonzeros
$N = 251001; #mat dim

$m1 = -1;
open MATVEC, "$matvec";
open PP,">$mat";
print PP $N, "\n";
$line = <MATVEC>;
$line = <MATVEC>;
$line = <MATVEC>;
$line = <MATVEC>;
for (my $k=0; $k<$n; $k++) { 
    $line = <MATVEC>;
    my @ligne = split " ",$line;
    if (@ligne>=0) {
	print PP $ligne[0]-1,' ',$ligne[1]-1, ' ',$ligne[2], "\n";
    }
}
print PP $m1,' ',$m1,' ',$m1,"\n";
close PP;

open PP,">$vec";
$line = <MATVEC>;
$line = <MATVEC>;
$line = <MATVEC>;
for (my $k=0; $k<$N; $k++) { 
    $line = <MATVEC>;
    my @ligne = split " ",$line;
    if (@ligne>=0) {
	print PP $ligne[0], "\n";
    }
}
close PP;
close MATVEC;
