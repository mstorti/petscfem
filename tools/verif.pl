#!/usr/bin/perl

$file = shift();
$t = shift();

sub add {
    my ($e,$t,$p,$l) = @_;
    $type{".$e"} = $t;
    $type{"->$t"} = $p;
}

__END__
add('m','octave','#.*\$Id: verif.pl,v 1.1 2002/02/25 23:09:46 mstorti Exp $Id: verif.pl,v 1.1 2002/02/25 23:09:46 mstorti Exp $ ');

if ($t != undef) {
    ($p,$l) = $type{'->'.$type};
    die "type $t not in table\n" unless $p;
    last;
} 
if (!$p && $file) {
    $file =~ /\.(\w*)$/;
    $e = $1;
#    if ($e && $t = $type{".$e"}) {
    if ($e) {
	($p,$l) = $type{$t};
    }
} 

die "couldn't find file: \"$file\", type: \"$type\"\n" unless $p;

open FILE,"$file";
while (<FILE>) {
    if ($p) {
	print "CVS entry OK: \"$_\"";
	last;
    }
}
