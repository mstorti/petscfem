#!/usr/bin/perl
# $Id: insdeb.pl,v 1.2 2002/01/14 03:45:06 mstorti Exp $

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

$ver = ' $Id: insdeb.pl,v 1.2 2002/01/14 03:45:06 mstorti Exp $ ';
$ver =~ s/\$/%/;
print OUT "// Insert Debug Info for SMC Finite State Machine generated code\n",
    "// script version: $ver\n";
    "// This key prevents multiple procesing: $key\n";

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
	/`/; print OUT <<EOF;
  if (s.matrix_p->print_fsm_transition_info_f())
    printf("from: \\"$from\\", event: \\"$event\\", "
           "to: \\"$to\\"\\n");
EOF
    /`/;
    flushin();
    } else {
	print OUT;
    }	
}
flushin();

close IN;
close OUT;
