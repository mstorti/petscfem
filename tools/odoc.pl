#!/usr/bin/perl
#__INSERT_LICENSE__
# $Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

@odoc=();

use Getopt::Std;
use English;
use Data::Dumper;

getopts("wWs:o:he:C:f:");

my $wiki_syntax = 1;
my $include_file = 1;
my $verb_char = '+';
$include_file = $opt_f if defined $opt_f;

$wiki_syntax = 0 if $opt_W;
$wiki_syntax = 1 if $opt_w;

#  if ($opt_h) {
#  /`/;
#    print <<'EOM';

sub current_file {
    $file = $ARGV;
    $file = $POSTMATCH if $opt_C && $ARGV =~ /^$opt_C/;
    return $file;
}

$tgetopt_pat='(?:T|E)GETOPT\w*\(\S*,(\w*),(\w*),([^ ,]*)\)';
$getopt_pat='GETOPT\w*\((\w*),(\w*),([^ ,]*)\)';

$otarget = ($opt_s ? $opt_s : "");

$sep="%".("----<>" x 10)."\n";
sub get_section {
    my $texfile = shift();
    my $section = shift();
    die "couldn't open $texfile\n" unless open TEX,$texfile;
    my @doc=();
    while (<TEX>) {
	goto FOUND if /^%section\s*$section\s*$/;
    }
    die "coudln't find section \"$section\" in file \"$texfile\"\n";
  FOUND:;
    while (<TEX>) {
	last  if /%section /;
	goto END_FOUND if /^%end_section\s*$/;
#	chomp;
	push @doc,$_;
    }
    die "last read: $_ Couldn't find \"%end_section\" tag\n";
  END_FOUND:;
    close TEX;
    return \@doc;
}

sub wiki2 {
    my ($ref,$wc,$pre,$post) = @_;
    my $text = $$ref;
    my $processed = "";
    while ($text =~ /(\s)$wc(\S)$wc(\s|[,.])/) {
	$processed .= "$`$1$pre$2$post";
	$text = "$3$'";
    }
    $text = "$processed$text";
    $processed = "";
    while ($text =~ /(\s)$wc(\S.*?\S)$wc(\s|[,.])/) {
	$processed .= "$`$1$pre$2$post";
	$text = "$3$'";
    }
    $$ref = "$processed$text";
}

sub wiki {
    my $t = shift();
    wiki2(\$t,"#","\\verb$verb_char","$verb_char");
    wiki2(\$t,"\\*","\\textbf{","}");
    wiki2(\$t,"_","\\emph{","}");
    return $t;
}

sub wiki_texi {
    my $t = shift();
    wiki2(\$t,"#","\@samp{","}");
    wiki2(\$t,"\\*","\@strong{","}");
    wiki2(\$t,"_","\@emph{","}");
    return $t;
}

$otargetf="";
@doclist=();
while (<>) {

    $otarget=$1 if m|//target (\s*)|;
    $otarget="" if m|//end_target|;
    $wiki_syntax = 1 if m|//__ENABLE_WIKI__|;
    $wiki_syntax = 0 if m|//__DISABLE_WIKI__|;
    $verb_char = '|' if m|//__USE_PIPE_FOR_VERB_CHAR__|;
    next if $otarget ne $otargetf;
    if ((/$tgetopt_pat/o || /$getopt_pat/o)
	&& ! m|//nd|) {
	chomp;
	print "in $ARGV line $.: line \"$_\" has no previous doc section\n"; 
    }
    if (m|^\s*//o (.*)|) {
	@doc=();
	push @doc,"$1\n";
	while (<>) {
	    if (m|^\s*// (.*)|) {
		push @doc,"$1\n";
		last if /_END\s*$/;
	    } elsif (m|^\s*//i_tex\s*(\S*)\s*(\w*)|) {
		my $texfile=$1;
		my $section=$2;
		# If not an absolute path then prepend the current
		# file directory
		if ($texfile !~ m|^/| && $ARGV =~ m|^(.*/)[^/]*|) {
		    $texfile = "$1$texfile";
		}
		my $tex = get_section($texfile,$section);
		my $l;
		push @doc,@{$tex};
	    } else {
		last;
	    }
	}

	if ($doc[0] =~ /\s*_T:(.*)/s) {
	    @docf=@doc;
	    $doc=join("",@doc);
	    die "in $ARGV. near line $.: no \"_T:\" tag in explicit doc: \n",
	    @docf unless $doc=~/_T:(.*\s)_N:/s;
	    $type=$1;
	    $doc=$POSTMATCH; 
	    $type =~s/\n/ /g;
	    $type =~ s/\s*$//;
	    die "in $ARGV near line $.: ",
	    "no \"_D:\" tag in explicit doc: \n",@docf,"\n-- doc here: >>  $doc",
	    @docf unless $doc=~/_D:/s;
	    $name=$`;
	    $doc=$POSTMATCH;
	    $name =~s/\n/ /g;
	    $name =~ s/\s*$//;
	    die "in $ARGV near line $.: ",
	    "no \"_DOC:\" tag in explicit doc: \n",@docf unless $doc=~/_DOC:/s;
	    $default=$`;
	    $doc=$POSTMATCH; 
	    $default =~s/\n/ /g;
	    $default =~ s/\s*$//;
	    die "in $ARGV near line $.: ",
	    "no \"_END\" tag in explicit doc: ",
	    @docf unless $doc=~/_END\s*$/s;
	    $doc=$`;
            $file = current_file();
            my $orig = ($include_file 
                        ? " (found in file: \\verb+$file+)\n"
                        : "");
	    $text = join("","\\item\\verb+$type+ \\verb+$name+ ",
			 "{\\rm(default=\\verb|$default|)}:\n",
			 $doc,$orig,$sep);
	    $name = $1 if $name =~ /^\s*(.*)\s*$/;
	    push @doclist,[$name,$type,$default,$text,$file,$doc,
			   $wiki_syntax];
	} elsif (/$tgetopt_pat/o || /$getopt_pat/o) {
	    $type=$1;
	    $name=$2;
	    $default=$3;
	    $tname=$name;
	    $tname =~ s/_/\\_/g;
	    $tdefault = $default;
	    $tdefault =~ s/!/!!/;
     my $header = <<EOT;
\\index{$tname@\\verb+$name+}
\\item\\verb+$type $name+ {\\rm(default=\\verb|$tdefault|)}:\n
EOT
	    $doc = join("",@doc);
            my $file = current_file();
            my $orig = ($include_file 
                        ? " (found in file: \\verb+$file+)\n"
                        : "");
	    $text = join("",$header,$doc,$orig,$sep);
	    $name = $1 if $name =~ /^\s*(.*)\s*$/;
	    push @doclist,[$name,$type,$tdefault,$text,$file,$doc,
			   $wiki_syntax];
#  	    print "doc: ",join("\n",@doc),"\n";
#  	    print "type: $type, name: $name, def: $default\n\n";
#  	} elsif (m|^\s*//e .*$|) {
#  	    eval $1;
#  	    print "doc: ",join("\n",@doc),"\n";
#  	    print "type: $type, name: $name, def: $default\n\n";
	} else {
	    print "GETOPT section not found after doc: ",join("",@doc),"\n";
	}
    }
}

@doclist = sort { lc($a->[0]) cmp lc($b->[0]) } @doclist;

#foreach my $r (@doclist) { print "$r->[0]\n"; }

if ($opt_o) {
    die "couldn't open \"$opt_o\"\n" unless open TEXOUT,">$opt_o";
$warn = <<EOT;
%-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>
% DONT EDIT MANUALLY THIS FILE !!!!!!
% This files automatically generated by odoc.pl from 
% source file \"$texfile\"
%-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>
EOT

    print TEXOUT $warn;
    foreach $doc (@doclist) {
	my $wiki = $doc->[6];
	my $t = $doc->[3];
	$t = wiki($t) if $wiki;
	print TEXOUT "$t";
    }
    print TEXOUT $warn;
    close TEXOUT;
}

## No more used now an info file is generated below
if (0 && $opt_e) {
    my %h = ();
    foreach $doc (@doclist) { $h{$doc->[0]} = 1; }
    my @u_docs = keys(%h);
    undef %h;
    @u_docs = sort @u_docs;
    die "couldn't open $opt_e\n" unless open EOUT,">$opt_e";
    print EOUT "(setq petscfem-option-list '(\n";
    foreach $doc (@u_docs) { print EOUT "(\"$doc\")\n"; }
    print EOUT "))\n";
    close EOUT;
}

if ($opt_e) {
    die "can't open \"$opt_e\"" unless open TEXI,">$opt_e";
    print scalar(@doclist),"\n";
    while (@doclist) { 
	my $k = $doclist[0]->[0];
	my @doc_text = ("\@node $k\n",
			"\@section $k\n",
			"\@vindex $k\n\n");
	while (@doclist && $doclist[0]->[0] eq $k) {
	    my $doc = shift @doclist;
	    my $wiki = $doc->[6];
	    my $d = $doc->[5];
	    $d =~ s/@/@@/g;
	    $d =~ s/\}/@\}/g;
	    $d =~ s/\{/@\{/g;
	    $d = wiki_texi($d) if $wiki;
	    my $node = $k;
	    push @doc_text, 
	    "$d\n",
	    "\@noindent\n",
	    "[Type: $doc->[1]]\@*\n",
	    "[Default: $doc->[2]]\@*\n",
	    "[Found in file: \"$doc->[4]\"]\@*\n",
	    "[Yank: <<$k>>]\n",
	    "\@c -----------------------------------\n";
	}
	print TEXI @doc_text;
    }
    close TEXI;
}

__END__

=head1 NAME 

odoc.pl - generate documentation for PETSc-FEM options in text-hashes

=head1 SYNOPSIS

    $ odoc.pl [OPTIONS] FILES ...

=head1 DESCRIPTION

Options  to programs,  elemsets and  other parts  of  C<PETSc-FEM> are
passed via text-hashes, which are  then queried by the code via calls
to functions in the C<getopt> package. This is normally done via calls
to  the *GETOPT*  macros.  You  can document  the  options so  defined
introducing special comments, for instance

  //o Use the weak form for the Galerkin part of the advective term. 
  SGETOPTDEF(int,weak_form,1);

will generate a LaTeX entry in the form 

    * int weak_form (default: 1) : Use the weak form for the Galerkin 
                       part of the advective term. 

You  can  put  almost any  kind  of  LaTeX  material in  the  embedded
documentation.  You must add a  doc section for each *GETOPT* call and
then execute

  $ odoc.pl -o doc.tex file1 file2 ...

in order to obtain a F<doc.tex> file. 

=head1 SYNTAX

The syntaxis of the embedded documentation is as follows:

=over
=item *

The  embedded  documentation  must  be  in a  block  B<preceding>  the
C<*GETOPT*> call.

=item * 

B<Spanning doc  blocks over more  than one line:> All  lines following
the  the C<//o>  directive  must start  with  C<// >  (i.e. the  C<//>
comment  directive must  be followed  by a  space.  Lines  starting by
C<//> followed by a non space character are 'special directives'.

=item *

B<The  >C<//i_tex>B<  directive  -  Including large  amount  of  LaTeX
material:> If  you have a too  large section then you  can introduce a
special line with a directive of the form C<//i_tex> as, for instance

  //o Sets the frequency save for the ``rotary save'' mechanism. 
  //i_tex nsdoc.tex rotary_save
  GETOPTDEF(int,nsaverot,100);

This  will read  the  section C<rotary_save>  entry  in the  nsdoc.tex
file. The  section in  the F<nsdoc.tex> file  are the lines  between a
line containing C<%section rotary_save> and C<%end_section>

   .... % other sections 

   %section rotary_save
   This section shows how to ...
   ...
   %end_section

   .... % other sections 

=item *

B<Including  explicit doc sections:>  If the  C<*GETOPT*> call  is too
complicated or the  C<getopt> functions (as C<get_int>, C<get_double>,
etc...) have been called explicitly, then you can include explicit doc
sections where you explicitly give the different fields of the entries
(C<type>,  C<name>,  C<default>,  C<doc>).  You  enter  them  in  the
following way

    //o _T: double[ndim]/double[ndim*ndof]/double[ndim*ndof*ndof] 
    //  _N: advective_jacobians _D: no default  _DOC: 
    //i_tex advdif.tex advective_jacobians
    //  _END

i.e.,  they are  delimited  by the  I<magic  strings> C<_T:>,  C<_N:>,
C<_D:>,  C<_DOC:>, C<_END>,  B<in that  order>.  You  can  put several
fields  in the  same  line.  You  can  use the  C<//i_tex> inside  the
C<_DOC:> block.

=item *

B<Not documented >C<*GETOPT*>B< calls:> If C<odoc.pl> finds a C<*GETOPT*>
call that doesn't have a previous C<//o> block, a warnng is issued. If
you want to suppress this  warning add a C<//nd> after the C<*GETOPT*>
call as, for instance

   TGETOPTDEF(GLOBAL_OPTIONS,double,alpha,1.); //nd

=item *

B<Wiki syntax:> A set of words between two stars (as in C<*foo bar
zoo*>) is bolded in the outoput (using C<\textbf{}>,
i.e. C<\textbf{foo bar zoo}>). Also, C<_foo bar_> expands to italic
(using C<\emph{}>) and C<#bar zoo#> expands to monospace (using
C<\verb++>). This is inspired in the syntax prevalent in most wiki
clones. (see for instance L<http://twiki.org>)

Several restrictions apply, however:

=over

=item * 

The whole construt must be contained in the same line, i.e. you can't
use them across lines.

=item * 

The character B<after> the first wiki character (i.e. one of C<*#_>)
and that one B<before> the second one must B<not> be white-space. The
characters before and after them B<must> be white-space. For instance
C<bar *foo zoo* mom> is expanded while C<bar * foo zoo * mom> and
C<bar*foo zoo* mom> not. The exact Perl regexp's are
C<\s$wc\S.*?\S$wc\s> and also C<\s$wc\S$wc\s> where C<$wc> is one of
the previously mentioned wiki characters.

=item * 

You can deactivate this feature using the special
C<//__DISABLE_WIKI__> and C<//__ENABLE_WIKI__> commands. 

=item * 

If the string to be expanded with C<#...#> construct contains a plus
sign C<+> you loose, because by default the string is expanded with 
the C<\verb+...+> construct. In that case you can include the 
command C<//__USE_PIPE_FOR_VERB_CHAR__>, and then the C<\verb|...|>
contruct will be used instead. 

=back

Example:

 //o This is an _a priori_ unstable method. *Do not* use
 //  unless you are sure what you do!! Use #soft_restart#
 //  instead. 

=back

=head1 OPTIONS

=over 4 

=item -s section_name

Process the specified section name in
the C<*.cpp> source file. Sections are delimited by C<//target
section_name> and C<//end_target> lines. For instance:

  // Contents of file `myprg.cpp'
  //target common_options
  ....
  //end_target

  //target specific_option
  ....
  //end_target

then will be able to write

   $ odoc.pl -s common_options myprg.cpp

=item -o outputfile

Put the generated LaTeX documentation in file C<outputfile>. 

=item -h

Give help.

=back

=head1 AUTHOR

Mario A. Storti E<lt>mario.storti@gmail.comE<gt>

=cut


