#!/usr/bin/perl

sub iso2ident {
    my @line = split "",shift();
    for (@line) {
	$_ = (/\w/ ? $_ : sprintf("_%02x",ord($_)));
    }
    return join("",@line);
}

sub ident2iso {
    my @line = split "",shift();
    my @out=();
    while (my $c= shift @line) {
	$c = chr(hex(shift(@line).shift(@line))) if $c eq '_';
	push @out,$c;
    }
    return join("",@out);
}

if ($0 eq 'iso2ident') {
    return iso2ident(shift());
} else {
    return ident2iso(shift());
}

=head1 NAME

iso2ident, ident2iso: encode/decode a string with alphanumeric chararacters 

=head1 SYNOPSIS

given a string C<iso2ident> converts it to a string composed of only
the C<[a-zA-Z0-9_]> set (i.e. alphanumeric + C<_>). Useful for
encoding arbitrary texts in filenames or identifiers.

=head1 DESCRIPTION
=head1 OPTIONS
=head1 RETURN VALUE
=head1 ERRORS
=head1 DIAGNOSTICS
=head1 EXAMPLES
=head1 ENVIRONMENT
=head1 FILES
=head1 CAVEATS
=head1 BUGS
=head1 RESTRICTIONS
=head1 NOTES
=head1 SEE ALSO
=head1 AUTHOR
=head1 HISTORY





perl

=cut  
# just for testing
while (<STDIN>) {
    chomp;
    my $t = iso2ident($_);
    print "encoded: \"$t\"\n";
    print "redeencoded: \"",ident2iso($t),"\"\n";
}
=cut
