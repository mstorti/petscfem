#!/usr/bin/perl
#__INSERT_LICENSE__

$pattern=shift();
$license_file=shift();
$file = shift();

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

@license=();
open LICENSE,$license_file;
@license=<LICENSE>;
close LICENSE;

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
