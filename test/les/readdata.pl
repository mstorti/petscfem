#!/usr/bin/perl

while (<>) {
    next unless /\[/;
    if (/^% Size \(\S*\) \(\S*\)/) {
	$n=$1; $m=$2;
	<>; <>; <>; 
	@mat=();
	while (<>) {
	    last if /\]/;
	    push @mat,$_;
	}
	$line = <>;
	die "no \"<ident> = spconvert(zzz);\" line found!\n"
	    unless $line =~ /^\s*\(\S*\)\s*= spconvert/;
	$name = $1;
	open MAT,">$name.mat.oct";
	print MAT @mat;
	close MAT;
    } else {
	
