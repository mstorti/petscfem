#__INSERT_LICENSE__
#$Id: myexpect.pl,v 1.9 2002/02/10 21:15:37 mstorti Exp $

# usage: expect($file,$pattern);
#Advances $file until finds a little of "\n" delimited $pattern's. 
$DEBUG_EXPECT = 0;

## position in count record
$OK= 1;
$NOT_OK= 2;
$CANT_OPEN= 3;

$COMPLAIN_ON_CANT_OPEN= 1 unless defined($COMPLAIN_ON_CANT_OPEN);

@stack=();
@output=();

sub read_file {
    my $file = shift();
    open FILE,$file;
    my $file = join "",(<FILE>);
    return $file;
}

sub printo {
    push @output,@_;
}

sub flush {
    print @output;
    @output=();
}

sub match_exactly {
    my ($line,$patt) = @_;
    $np = length($patt);
    if (length($line) < $np) {return 0;};
    return (substr($line,0,$np) eq $patt);
}

sub match_regexp {
    my ($line,$patt) = @_;
    return ($line =~ /$patt/);
}

sub expect {
    ($file,$descr,$pattern_list) = @_;
    open (SAL,$file) || do {printo "can't open $file\n" 
				unless ! $COMPLAIN_ON_CANT_OPEN;
			    cant_open(); return;};
    @sal=(<SAL>);
    close SAL;
    printo "Testing: \"$descr\" on file \"$file\"...";
    print "\n" if $DEBUG_EXPECT;
    @pattern = split("\n",$pattern_list);
    $record=0;
    my $inc=1;
    my $skip=1;
    my $match_fun = \&match_regexp;
    while ($pattern=shift @pattern) {
	if ($pattern =~ /^__REWIND__$/) {
	    $inc=1;
	    $record=0;
	    next;
	}
	if ($pattern =~ /^__BACKWARD__$/) {
	    $inc=-1;
	    next;
	}
	if ($pattern =~ /^__FORWARD__$/) {
	    $inc=+1;
	    next;
	}
	if ($pattern =~ /^__SKIP__$/) {
	    $skip=1;
	    next;
	}
	if ($pattern =~ /^__NO_SKIP__$/) {
	    $skip=0;
	    next;
	}
	if ($pattern =~ /^__EXACT_MATCH__$/) {
	    $match_fun = \&match_exactly;
	    next;
	}
	if ($pattern =~ /^__REGEXP_MATCH__$/) {
	    $match_fun = \&match_regexp;
	    next;
	}
	print "trying pattern: \"$pattern\"...  \n" if $DEBUG_EXPECT;
	while ($record<=$#sal) {
	    $_ = $sal[$record];
	    $record += $inc;
	    print "$record: $_" if $DEBUG_EXPECT;
	    do {chomp; 
		printo "        -> found: \"$_\"\n" if $DEBUG_EXPECT; 
		goto NEXT;} if &{$match_fun}($_,$pattern);
	    last unless $skip;
	}
	not_ok();
	printo "not OK.\n        --->  Couldn't find \"$pattern\"\n";
	return;
      NEXT:;
    }
    printo "OK. \n";
    ok();
}

sub inc {
    my $f = shift();
    if ($#stack<0) {
	push @stack,["--- Tests ---",0,0,0];
    }
    my $t = pop @stack;
    $t->[$f]++;
    push @stack,$t;
}

sub ok { inc($OK);}
sub not_ok { inc($NOT_OK);}
sub cant_open { inc($CANT_OPEN);}

sub begin_section {
    my $s=shift();
    push @stack,[$s,0,0,0];
}

sub end_section {
    my $t = pop @stack;
    my $total = $t->[1]+$t->[2]+$t->[3];
    my $total_open = $t->[$OK] + $t->[$NOT_OK];
    if ($total_open || $COMPLAIN_ON_CANT_OPEN) {
	print "--------\nStart section: \"$t->[0]\"\n";
	flush();
	print "Summary: \"$t->[0]\"",
	" -- OK: $t->[$OK].  Not OK: $t->[$NOT_OK]. ",
	"Couldn't open: $t->[$CANT_OPEN]. Total: $total\n";
    }
    if ($#stack>=0) {
	my $tt = pop @stack;
	$tt->[$OK] += $t->[$OK];
	$tt->[$NOT_OK] += $t->[$NOT_OK];
	$tt->[$CANT_OPEN] += $t->[$CANT_OPEN];
	push @stack,$tt;
    }
}

sub final_check {
    while ($#stack>=0) {
	end_section();
    }
}

=head1 NAME

myexpext.pl: verifies output from program tests

=head1 SYNOPSIS

C<expect(FILE,MESSAGE,patlist)> verifies that each pattern in
C<patlist> matches lines in FILE.

=head1 DESCRIPTION

The simplest way of using C<expect()> is to write a list of patterns
that should be found in the output of the test. Suppose the output
contains:

  Output run of program QBFG running on day Sun Apr  8 08:48:23 ART 2001
  ..... more lines here
  Total volume: 34.56
  more lines here...
  Total area: 23.43
  more lines here...
  Total impedance: 46.4
  more lines here...

You yount to chech that the figures in lines C<Total...> are precise
to the first digit. The ouput is in file C<QBFG.out> and you write a
small perl program like this

   #!/usr/bin/perl
   
   require 'myexpect.pl';
   
   expect("QBFG.out","Check on ouput of QBFG.out",<<'EOT');
   Total volume: 34.5
   Total area: 23.4
   Total impedance: 4
   EOT

   final_check();

In the default mode, C<expect()> takes the first pattern at a time
and starts scanning the file from the beginning, each lie at a time
until it finds a line that matches the pattern. Patterns are the usual
Perl patterns. So that remember to escape asterisks 'C<*>', question
marks 'C<?>', dots and others. You can leave the dot unescaped since
it matches itself, but the pattern is less strict (dot matches any
other character also).  Normally, when entering patterns with a I<here
in> document, as in the previous example, you protect the backslash
characters in the pattern list using quotes in the C<'EOT'>
terminator.

If the pattern is not found an error is reported and the test is
counted as a failure. If a line matching is found, C<expect()> takes
the folowing pattern and continue scanning the file from the line
following the previous match. If all the patterns are matched, then
the test is counted as a succes. If C<FILE> can't be opened, then this
is reported separately from error. The final count is reported by a
call to C<final_check()>.

You can alter this default behavior adding I<magic tokens> in
the pattern list. The magic tokens are

=over 4

=item __REWIND__

Rewind the file, i.e. scan for the next match starting from the
beginning of the file, rather than from the last match. This is useful
when you don't know exactly the order in which the lines will appear.
For instance file

   #------ contents of file test1.out
   line at the beginning
   ...
   other line 
   ...

matches the following call

   expect("test1.out","Check on ouput of test1.out",<<'EOT');
   other line 
   __REWIND__
   line at the beginning
   EOT
  
thanks to the presence of the C<__REWIND__> directive. 

=item __BACKWARD__

Scan the file backward for the next and subsequent patterns.

=item __FORWARD__

Cancel the C<__BACKWARD__> directive and continue scanning forward. 

=item __NO_SKIP__

In no-skip mode the following pattern is found to match with exactly
the following line, rather to scan the file from the following line
down. 

=item __SKIP__

Return to skip mode. 

=item __EXACT_MATCH__

Pattern has to match exactly the beginning of the line. (No special 
meaning for characters like C<.>, C<*>, etc...

=item __REGEXP_MATCH__

Go back to regexp match. 

=back

=head2 Sections

Sometimes it is useful to divide tests into sections. For this call
C<begin_section("section name")> before each section, and
C<end_section()> at the end.  All enclosed calls to C<expect()> are
assumed to be in the same logical section of tests and a summary is
reported for that section.

=head1 AUTHOR

Mario A. Storti <mstorti@intec.unl.edu.ar>

=cut

1;    
