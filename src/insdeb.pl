#!/usr/bin/perl
# $Id: insdeb.pl,v 1.1.2.1 2002/01/09 12:26:30 mstorti Exp $

sub flushin { print OUT @in; @in = (); }

$key = "==insdebinfo==";
$in = shift();
$out = shift();

open IN,"$in";
while(<IN>) {
    if (/$key/) { print "FSM debug info already inserted\n"; exit; }
}

close IN;

open IN,"$in";
open OUT,">$out";

print OUT "// Insert Debug Info for SMC Finite State Machine generated code\n",
    "// $key -- $0\n";

@in = ();
while(<IN>) {
    if (/^void.*FSM(.*)State::(\S*)\(/) {
	$from = $1;
	$event = $2;
	print OUT;
	$next = <IN>;
	push @in,$next;
	if ( $next =~ /\s*\{/ ) { flushin(); next; }
	die "couldn't read to-state in \"$next\"" 
	    unless $next =~ /\s*s\.SetState\(\S*::(.*)State\)/;
	$to = $1;
	print "Inserting debug info: from \"$from\" to \"$to\"\n";
	print OUT 
	    "  printf(\"from: \\\"$from\\\", event: \\\"$event\\\",\"\n",
	    "         \"to: \\\"$to\\\"\\n\");\n";
	flushin();
    } else {
	print OUT;
    }	
}
flushin();

close IN;
close OUT;
