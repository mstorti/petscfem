#!/usr/bin/perl
#__INSERT_LICENSE__

open OUT,"turbchan.out";
while (<OUT>) {
    ($res) = /time_step 10,.*=\s*(\S*)/;
    last if $res;
}
close OUT;

print "Residual at iteration 10 = ",$res,", < 2e-2 OK : ",$res<2e-2,"\n";

open STATE,"state.out.tmp";
$symm=1;
$v_small=1;
while (<STATE>) {
    ($u,$v,$h,$k,$e) = split " ",$_;
    <STATE>;
    ($u2,$v2,$h2,$k2,$e2) = split " ",$_;
    if (abs($v)>1e-10 || abs($v2)>1e-10) {$v_small = 0;}
    if (abs($u2-$u)>1e-10 || abs($h2-$h)>1e-10 ||
	abs($k2-$k)>1e-10 ||  abs($e2-$e)>1e-10) {$symm=0;}
}

print "symmetry OK : ",$symm,"\n";
print "transversal velocity small OK : ",$v_small,"\n";
print "k at output ",$k,". Correct value OK : ",abs($k-3.7e-2)<1e-3,"\n";
print "e at output ",$e,". Correct value OK : ",abs($e-1.5e-2)<1e-3,"\n";
