#!/usr/bin/perl
#__INSERT_LICENSE__

open SAL,"salida.sal";
open ID,">id.dat";
while (<SAL>) {
    $neq = $1 if /Total number of degrees of freedom neq:\s*(\d*)/;
    print ID "$1 $2\n" if /row (\d*): \((\d*) -> 1.000000\)/ && $2<=$neq;
}
close SAL;
close ID;

open OUT,"datos.dat";

sub read_mat {
    open MAT,">>".shift();
    while (<OUT>) {
	last if /\[/;
    }
    
    while (<OUT>) {
	last if /\]/;
	print MAT;
    }
    print MAT " -1 -1 0.\n";
}

unlink "mat.dat";
read_mat("mat.dat");
read_mat("mat.dat");
read_mat("mat.dat");
read_mat("mat.dat");
    
close OUT;
