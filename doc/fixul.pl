#!/usr/bin/perl -w
#__INSERT_LICENSE__
# $Id: fixul.pl,v 1.3 2003/02/10 14:54:54 mstorti Exp $

# It seems that Latex2html has a bug and it doesn't close the lists
# with the partial contents for the childs, when splitting a document. 
# This simple script fixes this, by counting the number of <UL> lines and
# of </UL> lines, and adding as many </UL> lines as needed in order
# to balance the tags. 

use strict;

my $name = shift();
my $open = 0;
my $found = 0;

open IN,"$name";
while(<IN>) {
    $open++ if m|^<UL>|i;
    $open-- if m|^</UL>|i;
    $found = 1 if m/^<LI> <A NAME=\"tex2html/;
}
close IN;
exit 0 if !$open;

print "$name: open $open, fixing...\n" if $open || $found;

my $bak = "$name.bak";
unlink $bak if -f $bak;
rename $name, $bak;
open IN,"$bak";
my @lines;
my ($last,$j);
$j = 0;
while (<IN>) {
    push @lines,$_;
    $last = $j if /^<LI>/i;
    $j++;
}
print "last <li> line $last: \"$lines[$last]\"\n";
close IN;
die "can't find <li> line\n" unless defined $last;
while ($open--) { $lines[$last+1] = "</UL>\n".$lines[$last+1]; }

open FIXED,">$name";
print FIXED @lines;
close FIXED;
