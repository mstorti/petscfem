#!/usr/bin/perl

use Getopt::Std;

getopts("hsm:");

$hdrmxlen = 22;
$sep = " " x $hdrmxlen;
$MGS=($opt_m ? $opt_m : ""); # magic prefix
$hdrpat = "'\\\$ make \$target' : ";

sub printdoc {
    # Sort 
    if ($opt_s) {
	@docl = sort {$a->[0] cmp $b->[0]} @docl;
    }
    # Print entries
    while ($doci = shift @docl) {
	@doc = @{$doci};
	$target= shift @doc;
	eval "\$hdr = \"$hdrpat\"";
	$nh = length($hdr);
	print $hdr;
	if ($nh<$hdrmxlen) {
	    print " " x ($hdrmxlen-$nh);
	}
	print shift(@doc);
	foreach $l (@doc) { 
	    $l =~ s/^\s*//;
	    $l =~ s/^>//;
	    print $sep,$l;
	}
    }
}

if ($opt_h) {system "pod2text $0"; exit;}

@doclist=();

@doc = ();
$target="";

while (<>) {
    if (/^\#${MGS}s (.*)$/o || /^\#${MGS}s(\s*)$/o) {
	printdoc if $#docl>=0;
	print "-" x 80,"\n  ----------- $1 ------------\n",
	"-" x 80,"\n" unless $1=~/^\s*$/;
        @doc=();
    } elsif (/^\#${MGS}w (.*)$/o) {
	push @doc,"$1\n";
    } elsif (/^\#${MGS}v (.*)$/o) {
        eval $1;
    } elsif (/^(\S*)\s*:/) {
	next unless $#doc>=0;
	$target = $1;
	push @docl,[$target,@doc];
	@doc=();
    } elsif (/^\#${MGS}e (\S*)\s*$/o) {
	$target = $1;
	push @docl,[$target,@doc];
	@doc=();
    } elsif (/^\#${MGS}p$/o) {
	print;
    } elsif (/^\#${MGS}p (.*)$/o) {
	print "$1\n";
    } else {
	@doc = ();
    }
}

printdoc if $#docl>=0;

=head1 NAME

    makewhat.pl - Prints info on targets extracted from Makefiles

=head1 SYNOPSIS

  $ makewhat.pl [options] <makefiles> ... 

=head1 OPTIONS

=over 4

=item -h

give help

=item -s

Sort entries alphabetically in a given section

=item -m <magic_prefix>

Set the magic prefix. 

=back

=head1 DESCRIPTION

Add a comment just above the corresponding target in a Makefile in the form

    #w Builds the library and 
    #w cleans the directory
    buildclean: 
            commands
            ...

Then when running C<makewhat.pl> it scans the makefiles, strips this
comments and prints somethinng like

    '$ make buildclean'   :  Builds the library and
			     cleans the directory
    '$ make other_target' :  Other description here ...

Usually, one inserts a target of the form

     #w Prints info on targets
     what:
            @makewhat.pl Makefile Makefile.base 

So that a user can do

    $ make what

    'make buildclean'  :  Builds the library a
			  cleans the directory

    'make other_target' :  Other description here ...

    'make what'  : Prints info on targets

Leading whitespace in C<#w> lines is removed. If you want to preserve
it, start with a 'E<gt>' sign. Leading whitespace until the 'E<gt>' is
removed and the line is taken literally beginning at the character
following the 'E<gt>' sign.

=head2 Other directives

=over

=item Sections

Lines of the form

  #s Section name

separate targets in 'sections'. If sorting is enabled (option C<-s>)
then targets are sorted within each section. A target with no section
name flushes the target list.

=item Textual lines

Lines starting with C<#p> are printed textually (the C<#p> is
stripped) and printed out as the file is parsed. This is used to print
banners and legends like C<[In file path/to/subdir/Makefile]>. 

=item Change target and nonexistent targets

A line starting with C<#e name> is equivalent to finding a target of
the form C<name:> so that a target message can be generated with a
block of comments of the form

   #w This is a dummy target
   #w and will be reported even if there is not
   #w a real target
   #e dummy

This would print

  '$ make dummy' : This is a dummy target
                   and will be reported even if there is not
                   a real target

This can be also used for replacing the name of the target. For instance targets of the form 

 #w Converts file.a to file.b
 %.b: %.a
         myconvert $< $@

would generate a message of the form 

 '$ make %.b' :        Converts file.a to file.b

It is more readable to use

 #w Converts file.a to file.b
 #e file.b
 %.b: %.a
         myconvert $< $@

which would generate

 '$ make file.b' :     Converts file.a to file.b

instead.

=item Insert Perl code

A line of the form C<#v perl-code> causes C<perl-code> to be evaled,
and allows the user to change some parmeters. For instance

=over

=item C<$hdrmxlen> [default 22]

The length of the target part of the message. Subsequent lines are indented this amount. 

=item C<$hdrpat> [default C<'\\\$ make \$target' : >]

The target part is printed with this pattern. Note that this string
has to be double escaped, since it is interpolated twice. 

=back

=back

=head2 Changing the magic prefix

If you think that the C<#letter> commands are too
loose, then you can add some magic prefix so that the commands are now
C<#{magic-prefix}letter>. For instance if you set C<-m MKW> in the
options then the commands are now C<#MKWw>, C<#MKWs>, etc... Note that
this string will be evaluated in a pattern environment so that escape
characters like C<*.?/> etc...

=head2 

=head1 AUTHOR

Mario A. Storti E<lt>mstorti@intec.unl.edu.arE<gt>

=cut
    
