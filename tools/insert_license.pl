#!/usr/bin/perl
#__INSERT_LICENSE__

use Getopt::Std;
getopts("u");

$license_file=shift();
$file = shift();

@license=();
open LICENSE,$license_file;
$pattern = <LICENSE>;
chomp $pattern;
$pattern =~ s/<L>/INSERT_LICENSE/g;
$tag  = <LICENSE>;
$tag =~ s/<L>/INSERT_LICENSE/g;
@license=<LICENSE>;
close LICENSE;

if (!$opt_u) {

    open INFILE,"$file";
    $matches=0;
    while (<INFILE>) {
	if (m|$pattern|) {
	    $matches=1;
	    last;
	}
    }
    close INFILE;

    exit(0) if !$matches;

    rename $file,"$file~";
    open INFILE,"$file~";
    open OUTFILE,">$file";
    while (<INFILE>) {
	if (m|$pattern| && !/__NO_INSERT__/) {
	    print OUTFILE @license;
	} else {
	    print OUTFILE $_;
	}
    }
    close INFILE;
    close OUTFILE;

    @l=stat("$file~");
    chmod $l[2],$file;

    unlink "$file~";

} else {

    open INFILE,"$file";
    @outfile = ();
    @buf=();
    $l = 0;
    $modif=0;
    while (<INFILE>) {
	if ($_ eq $license[$l]) {
	    $l++;
	    push @buf;
	} else {
	    $l=0;
	    push @outfile,@buf,$_;
	}
	if ($l == @license) {
	    $modif=1;
	    push @outfile,$tag;
	    @buf=();
	    $l=0;
	}
    }
    close INFILE;
    if ($modif) {
	open OUTFILE,">$file";
	print OUTFILE @outfile;
	close OUTFILE;
    }
}
