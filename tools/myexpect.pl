#__INSERT_LICENSE__
#$Id: myexpect.pl,v 1.19 2003/11/16 19:56:30 mstorti Exp $

use English;
## position in count record
$OK= 1;
$NOT_OK= 2;
$CANT_OPEN= 3;

$COMPLAIN_ON_CANT_OPEN= 1 unless defined($COMPLAIN_ON_CANT_OPEN);

@stack=();
@output=();
@stack_output=();

$PRE = "{{";
$POST = "}}";
$SEP = "}{";
$COMMENT = "#>>";
$WD = "";

sub P { print @_; }

sub read_file {
    my $file = shift();
    die "can't open file \"$file\"\n" unless  open FILE,$file;
    my $file = join "",(<FILE>);
    return $file;
}

sub printo_push {
    push @stack_output,[@output];
    @output = ();
}

sub printo_pop {
    my $last = pop @stack_output;
    flush() unless defined $last;
    unshift @output,@$last;
}

sub printo_discard {
    my $last = pop @stack_output;
    flush() unless defined $last;
    @output = ();
    unshift @output,@$last;
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
    ## Parse pattern
    my $orig_patt = $patt;
    my $new_patt = "";
    my @conds = ();
    while ($patt =~ /$PRE/) {
	$patt = $POSTMATCH;
	$new_patt .= $PREMATCH;
	die "delimiters doesn't match in pattern: \"$orig_patt\"\n" 
	    unless $patt =~/$POST/;
	$cond = $PREMATCH;
	$patt = $POSTMATCH;
	my $subpat = ".*";
	if ($cond =~ /$SEP/) {
	    # my $condo = $cond;
	    $cond = $POSTMATCH;
	    $subpat = $PREMATCH;
	    # print "condo: $condo, subpat: $subpat, cond: $cond\n";
	}
	$new_patt .= "($subpat)";
	push @conds,$cond;
    }
    $new_patt .= $patt;
    # print "new_patt: $new_patt\n";
    # print "conds: <",join("><",@conds),">\n";
    my @matches = ($line =~ /$new_patt/);
    ## return 0 unless the line matched
    return 0 unless @matches;
    ## return 0 if there are conditions and they didn't match
    ## in number with the matches
    return 0 if @conds && @matches != @conds;
    ## check each condition if any
    for (my $j=0; $j<@conds; $j++) {
	$W = $matches[$j];
	return 0 unless eval $conds[$j];
    }
    return 1;
}

sub read_output {
    my ($sal,$file) = @_;
    open (SAL,$file) || do {printo "can't open $file\n" 
				unless ! $COMPLAIN_ON_CANT_OPEN;
			    cant_open(); return 1; };
    @$sal=(<SAL>);
    close SAL;
    return 0;
}

sub expect {
    ($file,$descr,$pattern_list) = @_;
    @sal = ();
    return if $file && $file!~ m|/$| && read_output(\@sal,$file);
    $WD = $1 if $file =~ m|^(.*/)[^/]*$|;
    printo_push();
    printo "Testing: \"$descr\" on file \"$file\"...";
    printo "\n" if $opt_d;
    @pattern = split("\n",$pattern_list);
    $record=0;
    my $inc=1;
    my $skip=1;
    my $match_fun = \&match_regexp;
    my $file_changed = 0;
    while ($pattern=shift @pattern) {
	if ($pattern =~ /^$COMMENT/) { next; }
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
	if ($pattern =~ /^__CONFIG__/) {
	    $rest = $POSTMATCH;
	    while ($rest =~ /^\s*(\S*)\s*(\S)/) {
		my $var = $1;
		my $delim = $2;
		$rest = $POSTMATCH;
		die "not valid var in config line var \"$var\", line: \"$pattern\"\n"
		    unless $var =~ /^(PRE|SEP|POST|WD|COMMENT)$/;
		die "unbalanced delimiters at config line: \"$pattern\"\n"
		    unless $rest =~ /^(.*?)$delim/;
		my $value = $1;
		$rest = $POSTMATCH;
		${$var} = $value;
		## print "setting \$${var} = $value\n";
	    }
	    die "can't parse config line: \"$pattern\"\n" unless $rest =~ /^\s*$/;
	    next;
	}
	if ($pattern =~ /^__SWITCH_FILE__\s*(\S*)\s*$/) {
	    $file = $1;
	    $file = "$WD/$file" unless $file=~ m|^/|;
	    @sal = ();
	    read_output(\@sal,$file);
	    $record=0;
	    $file_changed = 1;
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
	printo "trying pattern: \"$pattern\"...  \n" if $opt_d;
	while ($record<=$#sal) {
	    $_ = $sal[$record];
	    $record += $inc;
	    printo "$record: $_" if $opt_d;
	    do {chomp; 
		printo "        -> found: \"$_\"\n" if $opt_d; 
		goto NEXT;} if &{$match_fun}($_,$pattern);
	    last unless $skip;
	}
	not_ok();
	printo "not OK.\n        --->  Couldn't find \"$pattern\"\n";
	printo "    [on file $file]\n" if $file_changed;
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

sub ok { 
    if ($opt_n) { 
	printo_discard(); 
    } else {
	printo_pop(); 
    }
    inc($OK);
}
sub not_ok { printo_pop(); inc($NOT_OK);}
sub cant_open { printo_pop(); inc($CANT_OPEN);}

sub begin_section {
    my $s=shift();
    push @stack,[$s,0,0,0];
}

sub end_section {
    my $t = pop @stack;
    my $total = $t->[1]+$t->[2]+$t->[3];
    my $total_open = $t->[$OK] + $t->[$NOT_OK];
    my $print_sec = !$opt_n && ($total_open || $COMPLAIN_ON_CANT_OPEN)
	|| $opt_n && $t->[$NOT_OK];
    if ($print_sec) {
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

You want to check that figures in lines C<Total...> are precise
to the first digit. The output is in file C<QBFG.out> and you write a
small perl program like this

   #!/usr/bin/perl
   
   require 'myexpect.pl';
   
   expect("QBFG.out","Check on output of QBFG.out",<<'EOT');
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
the following pattern and continue scanning the file from the line
following the previous match. If all the patterns are matched, then
the test is counted as a success. If C<FILE> can't be opened, then this
is reported separately from error. The final count is reported by a
call to C<final_check()>.

=head2 Special directives

You can alter this default behavior adding I<special directives> in
the pattern list. They are

=over 4

=item C<__REWIND__>

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

   expect("test1.out","Check on output of test1.out",<<'EOT');
   other line 
   __REWIND__
   line at the beginning
   EOT
  
thanks to the presence of the L<__REWIND__> directive. 

=item C<__BACKWARD__>

Scan the file backward for the next and subsequent patterns.

=item C<__FORWARD__>

Cancel the L<__BACKWARD__> directive and continue scanning forward. 

=item C<__NO_SKIP__>

In no-skip mode the following pattern line is checked to match with
exactly the following output line, rather to scan the whole output
file from the following line down.

=item C<__SKIP__>

Return to skip mode. 

=item C<__EXACT_MATCH__>

Pattern has to match exactly the beginning of the line. (No special 
meaning for characters like C<.>, C<*>, etc...

=item C<__REGEXP_MATCH__>

Go back to regexp match. 

=item C<__SWITCH_FILE__ <file>>

Stop reading this output file and switch over to C<file>. File 
name is taken absolute if starts with C</> otherwise relative to the 
current working directory. The working directory is taken from the
first file entered or changed with the C<<< __CONFIG_ WD "<dir>" >>> 
(directive. Example:
 
   __SWITCH_FILE__ foodir/foofile.out

=item C<__CONFIG__>

Allows changing configuration variables. Syntax is

      __CONFIG__  var "value" ...

Set configuration variable C<$var> to C<value>. Current configuration
variables are

=over 8

=item C<PRE>, C<SEP>, C<POST>

The patterns that matches for the start, end and separator 
of an embedded block. (See below)

=item C<COMMENT>

The pattern that, when at the start of a line means a
comment. Default: C<<< #>> >>>.

=item C<WD>

The working directory.

=back

The double quote delimiter around C<value> may be changed by any
non-blank character. Examples:

   __CONFIG__ PRE "<<" SEP "><" POST ">>" WD "foodir/bardir" COMMENT "#"
   __CONFIG__ PRE COMMENT |#:|

The last one sets comment start to C<#:>. 

=back

=head2 Embedded blocks

You can embed Perl code blocks in the pattern for checking the result
of the matches themselves. The syntax is C<{{CONDITION}}> or
C<{{PATTERN}{CONDITION}}>. In the first form the C<{{...}}> block is
replaced with C<(.*)>. If the output line matches, then every
C<CONDITION> is evaluated and the line matches if every condition
returns true. Inside the condition the special variable C<$W> takes 
the value of the matching string. Example: the following output lines

  Total mass: 0.356, 
    max density: 0.48956, min density: 0.001234

match the following pattern lines

  Total mass: {{abs($W-0.355)<2e-3}}, 
    max density: {{abs($W-0.4)<0.1}}, min density: {{$W<2e-3}}

In the second form the C<PATTERN> section allows to specify the pattern 
which replaces the block. Example: the output line

  Processed case: ABC890

matches

  Processed case: {{[A-Z]*}{$W eq 'ABC'}}{{\d*}{$W>800 && $W < 900}}

The syntax of the block may be changed with the L<C<__CONFIG__>>
directive. The corresponding variables are C<PRE>, C<SEP> y <POST>
(mnemonic: prefix, separator and postfix). Possible choices are

  #>> Following three examples use matching delimiters (like <>, () or[])
  #>> warning: angles (<>) may collide with comparison expressions
  __CONFIG__ PRE "((" SEP ")(" POST "))"
  __CONFIG__ PRE "<<" SEP "><" POST ">>"
  __CONFIG__ PRE "[[" SEP "][" POST "]]"
  #>> This is very simple
  __CONFIG__ PRE "{" SEP "," POST "}"
  #>> Another simple one
  __CONFIG__ PRE "<" SEP "|" POST "}"

In order to avoid collision you can increase the delimiter levels, e.g.

  #>> Very paranoid
  __CONFIG__ PRE "{{{{" SEP "}}{{" POST "}}}}"
  #>> Combines with colon
  __CONFIG__ PRE "<<:" SEP ":><:" POST ":>>"

=head2 Sections

Sometimes it is useful to divide tests into sections. Start sections
with C<begin_section(<section name>)>, and end with C<end_section()>.
All enclosed calls to C<expect()> are assumed to be in the same
logical section of tests and a summary is reported for that section.
Example:

 begin_section("Navier-Stokes tests");
 expect("NS/output1.txt","NS Test 1","NS Test 1 OK");
 expect("NS/output1.txt","NS Test 2","NS Test 2 OK");
 expect("NS/output1.txt","NS Test 3","NS Test 3 OK");
 end_section();

 begin_section("Electro-magnetic tests");
 expect("EM/output1.txt","EM Test 1","EM Test 1 OK");
 expect("EM/output1.txt","EM Test 2","EM Test 2 OK");
 expect("EM/output1.txt","EM Test 3","EM Test 3 OK");
 expect("EM/output1.txt","EM Test 4","EM Test 4 OK");
 end_section();

=head1 AUTHOR

Mario A. Storti <mstorti@intec.unl.edu.ar>

=cut

1;    
