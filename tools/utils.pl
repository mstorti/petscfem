# -*- perl -*-
#$Id: utils.pl,v 1.2 2004/08/23 02:25:50 mstorti Exp $

# From the Perl FAQ 4
# usage: $num = getnum(
sub getnum {
    use POSIX qw(strtod);
    my $str = shift();
    return unless defined $str;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    $! = 0;
    my($num, $unparsed) = strtod($str);
    if (($str eq '') || ($unparsed != 0) || $!) {
	return undef;
    } else {
	return $num;
    }
}

sub is_numeric { defined getnum($_[0]) }

# usage: $array2 = aload("myfile")
#   returns an array ref to a 2D array 
#  after, gts elements with $array->[$i][$j]
sub aload {
    my $file = shift();
    die "can't open $file\n" unless open F,$file;
    my @data;
    while(<F>) { 
	my @v = map { $_ = getnum($_); } split(" ",$_);
	push @data,[@v]; 
    }
    return \@data;
}

sub count_lines {
    my $file = shift();
    die "passed undefined value!!\n" unless defined $file;
    open FILE,$file;
    my $lines=0;
    while (<FILE>) { $lines++; }
    return $lines;
}

sub connect_to0 {
    my $file = shift();
    my $file_out = shift();
    if (!defined $file_out) {
	if (-f "$file~") { unlink "$file~"; }
	rename $file,"$file~";
	$file_out = $file;
	$file = "$file~";
    }
    open FILE,$file;
    open FILEO,">$file_out";
    my $lines=0;
    while (<FILE>) { 
	my @l = split " ";
	for $f (@l) { print FILE$f-1," "; }
	print FILEO "\n";
	$lines++;
    }
    close FILE;
    close FILEO;
}

## usage: insertbef <text>, <file>, <regexp>;
## Insert <text> in <file> before line containing <regexp>
sub insertbef {
    my ($text,$file,$pos_rgxp) = @_;
    rename $file, "$file~";
    open FOLD,"$file~";
    open FNEW,">$file";
    while (<FOLD>) {
	print FNEW "$text\n" if /^$pos_rgxp/;
	print FNEW;
    }
    close FNEW;
    close FOLD;
}

sub ask2 {
    my $ans;
    my $text = shift();
    $text = "$text (y/n [default y]) > ";
    while (1) {
	print $text; 
	$ans = <STDIN>;
	chomp $ans;
	if ($ans =~ /^\s*y\s*$/i) { return 1; }
	elsif ($ans =~ /^\s*n\s*$/i) { return 0; }
	else { print "incorrect answer: \"$ans\"\n"; }
    }
}

sub ask {
    my $ans;
    my $text = shift();
    print "$text (y/n [default y]) > ";
    $ans = <STDIN>;
    chomp $ans;
    die "Aborting... \n" if $ans =~ /n/i;
}

## Greps a file for a pattern
## usage grepf(<pattern>,<file>);
sub grepf {
    my ($pattern,$file) = @_;
    # print "pattern \"$pattern\", file \"$file\"\n";
    die "can't find $file\n" unless -f $file;
    open FILE,$file || die "can't open $file\n";
    my $found;
    while (<FILE>) {
	if (/$pattern/) {
	    $found = $_;
	    last;
	}
    }
    close FILE;
    return $found;
}

sub sysdbg() { print @_,"\n"; }

## Count how many occurrences of $pat are in $s. 
## usage: $count = count_match($s,$pat);
sub count_match {
    my ($s,$pat) = @_;
    my $rest = $s;
    my $count = 0;
    while ($rest =~ /$pat/) {
	$rest = $POSTMATCH;
	$count++;
    }
    return $count;
}

## usage: transformer($text,$matcher)
## Process TEXT applying the MATCHER function until all
## the string is processed. 
## $matcher should be a reference tu function that takes
## as argument a string and returns a boolean, the already processed
## part of the string, and the part to be processed yet.
## 
## For instance $matcher = sub {
##    $t = shift();
##    return 1,$PREMATCH.$subst,$POSMATCH if $text =~ /$patt/;
## }
##
## makes transformer($text,$matcher) to behave as s/$text/$patt/ 
sub transformer {
    my ($text,$matcher)  = @_;
    my @newtext = ();
    while (1) {
	my ($pre,$post) = &$matcher($text);
	last unless defined $pre;
	push @newtext,$pre;
	$text = $post;
    }
    push @newtext,$text;
    return join("",@newtext);
}

1;
