#!/usr/bin/perl

$in_data=0;
$name="";
$basename="tmp_petsc_data";
@namelist=();

while(<>) {
    if (/^(\S*) = \[/) {
	$in_data=1;
	$name=$1;
	open TMP,">$basename._";
    } elsif (/\]/) {
	$in_data=0;
	close TMP;
	rename "$basename._","$basename.$name";
	push @namelist,$name;
    } elsif ($in_data) {
	print TMP;
    } elsif (/\s*(\S*)\s*=\s*spconvert/) {
	$nameo=$name;
	$name=$1;
	pop @namelist;
	push @namelist,$name;
	rename "$basename.$nameo","$basename.$name";
    } elsif (/^%/) {
    } elsif (/^zzz = zeros/) {
    } else {
	print "????? -> <$_>\n";
    }
}
	
open TMP,">${basename}_script.m";
foreach $name (@namelist) {
    print TMP "load $basename.$name; $name=$basename; ",
    "unlink(\"$basename.$name\"); clear $basename;\n";
}
close TMP;
