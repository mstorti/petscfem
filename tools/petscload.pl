#!/usr/bin/perl
#__INSERT_LICENSE__

$in_data=0;
$name="";
$basename="tmp_petsc_data";
@namelist=();
%names=();

sub check_name {
#    print "checking $name\n";
    if ( $names {$name} ) {
#  	my $nameo = $name;
#  	if ($name =~ /^(.*)_(\d*)$/) {
#  	    $name = "$1_".($2+1);
#  	} else {
#  	    $name = "${name}_0";
#  	}
	my $new_name,$c=0;
	while (1) {
	    $new_name = "${name}_$c";
#	    print "trying $new_name\n";
	    last unless $names{$new_name};
	    $c++;
	}
	print "repeated name, replacing $name -> $new_name\n";
	$name = $new_name;
    }
    $names{$name}=1;
}

while(<>) {
    if (/^(\S*) = \[/) {
	$in_data=1;
	$name=$1;
	check_name();
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
	$names{$nameo} = 0;
	$name=$1;
	check_name();
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
