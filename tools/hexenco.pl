#!/usr/bin/perl

$ESCAPE_CHAR = '_';

sub iso2ident {
    my @line = split "",shift();
    for (@line) {
	$_ = ( /[a-zA-Z0-9]/ ? $_ : $ESCAPE_CHAR.sprintf("%02x",ord($_)));
    }
    return join("",@line);
}

sub ident2iso {
    my @line = split "",shift();
    my @out=();
    while (1) {
	$c = shift @line;
	last unless defined $c;
	$c = chr(hex(shift(@line).shift(@line))) if $c eq $ESCAPE_CHAR;
	push @out,$c;
    }
    return join("",@out);
}

if ($0 =~ /iso2ident/) {
    print iso2ident(shift());
} else {
    print ident2iso(shift());
}

=head1 NAME

iso2ident, ident2iso: encode/decode a string with alphanumeric chararacters 

=head1 SYNOPSIS

Given a string C<iso2ident> converts it to a string composed of only
the C<[a-zA-Z0-9_]> set (i.e. alphanumeric + C<_>). Useful for
encoding arbitrary texts in filenames or identifiers. C<ident2iso>
does the inverse task.

=head1 DESCRIPTION

Suppose you want to have files with names "case: a=0.2345; c=3.56;
b=4.89e-15". You can do this in Unix but spaces and other characters
may lead to problems. The problem is harder if you want identifiers
for, say, C<C>, C<Perl> or a similar language. This utility converts
all characters different from the set C<[a-zA-Z0-9]> to C<_+hexcode>,
so that you end up with only characters in the set C<[a-zA-Z0-9_]>.
Underscore acts as a 'escape character'.  Of course, the escape
character has to be also encoded. The escape character may be changed
by modifying the value of C<$ESCAPE_CHAR>.

=head1 EXAMPLES

 $ iso2ident 'case: a=0.2345; c=3.56; b=4.89e-15'
 case_3a_20a_3d0_2e2345_3b_20c_3d3_2e56_3b_20b_3d4_2e89e_2d15
 $ ident2iso case_3a_20a_3d0_2e2345_3b_20c_3d3_2e56_3b_20b_3d4_2e89e_2d15
 case: a=0.2345; c=3.56; b=4.89e-15
 $ 

=head1 AUTHOR

Mario A. Storti C<mstorti@intec.unl.edu.ar>

=cut

=cut  
# just for testing
while (<STDIN>) {
    chomp;
    my $t = iso2ident($_);
    print "encoded: \"$t\"\n";
    print "redeencoded: \"",ident2iso($t),"\"\n";
}
=cut
