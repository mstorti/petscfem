#!/usr/bin/perl -w
use strict;

my $file = shift();
exit 0 if $file=~ /~$/;
my $is_func = 0;
my $has_license=0;

open FILE,$file;
while (<FILE>) { exit 0 if /__INSERT_LICENSE__/ || /^\s*function\s/; }
close FILE;

sub license {
    print NEWFILE <<'EOF';
##__INSERT_LICENSE__
## $Id: addheader.pl,v 1.1 2003/01/08 15:54:26 mstorti Exp $
EOF
}

rename $file,"$file~";
open FILE,"$file~";
open NEWFILE,">$file";
my $line = 0;
my $skip_first_comment_state=0;
while (<FILE>) {
    if (!$is_func && $line==0) { license(); }
    print NEWFILE;
    $line++;
}
close FILE;
close NEWFILE;
