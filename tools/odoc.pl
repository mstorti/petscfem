#!/usr/bin/perl

@odoc=();

use Getopt::Std;
getopts("s:o:h");

#  if ($opt_h) {
#  /`/;
#    print <<'EOM';

$tgetopt_pat='TGETOPT\w*\(\S*,(\w*),(\w*),([^ ,]*)\)';
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

$otargetf="";
@doclist=();
while (<>) {
    $otarget=$1 if m|//target (\s*)|;
    $otarget="" if m|//end_target|;
    next if $otarget ne $otargetf;
    if ((/$tgetopt_pat/o || /$getopt_pat/o)
	&& ! m|//nd|) {
	chomp;
	print "line \"$_\" has no previous doc section\n"; 
    }
    if (m|^\s*//o (.*)|) {
	@doc=();
	push @doc,$1;
	while (<>) {
	    if (m|^\s*// (.*)|) {
		push @doc,$1;
		last if /_END\s*$/;
	    } elsif (m|^\s*//i_tex\s*(\S*)\s*(\w*)|) {
		$texfile=$1;
		$section=$2;
		$tex = get_section($texfile,$section);
		push @doc,@{$tex};
	    } else {
		last;
	    }
	}
#  	if (/GETOPT.*_ND\(\S*,(\w*),(\w*),(\S*)\)/ 
#  	    || /GETOPT.*\((\w*),(\w*),(\S*)\)/) {
	if ($doc[0] =~ /\s*_T:(.*)/) {
	    @docf=@doc;
	    $doc=join("",@doc);
	    die "not \"_T:\" tag in explicit doc: \n",
	    @docf unless $doc=~/_T:(.*) _N:/;
	    $type=$1;
	    $doc=$';
	    die "not \"_D:\" tag in explicit doc: \n",@docf unless $doc=~/(.*) _D:/;
	    $name=$1;
	    $doc=$';
	    die "not \"_DOC:\" tag in explicit doc: ",@docf unless $doc=~/(.*) _DOC:/;
	    $default=$1;
	    $doc=$';
#	    print "doc antes de quitar text: $doc\n";
	    die "not \"_END\" tag in explicit doc: ",
	    @docf unless $doc=~/_END\s*$/;
	    $doc=$`;
#	    print "text in doc: $doc\n";
	    $text = join("","\\item\\verb+$type $name+ ",
			 "{\\rm(default=\\verb|$default|)}:\n",
			 $doc,$sep);
	    push @doclist,[$name,$type,$default,$text];
	} elsif (/$tgetopt_pat/o || /$getopt_pat/o) {
	    $type=$1;
	    $name=$2;
	    $default=$3;
	    $tname=$name;
	    $tname =~ s/_/\\_/g;
	    /`/;
	    $text = join("",<<EOT,join("\n",@doc),"\n",$sep);
\\index{$tname@\\verb+$name+}
\\item\\verb+$type $name+ {\\rm(default=\\verb|$default|)}:\n
EOT
/`/;
	    push @doclist,[$name,$type,$default,$text];
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

@doclist = sort {$a[0] cmp $b[0]} @doclist;

die "couldn't open $opt_o\n" unless open TEXOUT,">$opt_o";
/`/;
$warn = <<EOT;
%-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>
% DON'T EDIT MANUALLY THIS FILE !!!!!!
% This files automatically generated by odoc.pl from 
% source file \"$texfile\"
%-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>
EOT
/`/;

print TEXOUT $warn;
foreach $doc (@doclist) {
    print TEXOUT $doc->[3];
}
print TEXOUT $warn;
close TEXOUT;

__END__

=head1 NAME 

odoc.pl - generate documentation for PETSc-FEM options in text-hashes

=head1 SYNOPSIS

    $ odoc.pl [OPTIONS] FILES ...

=head1 DESCRIPTION

Options  to programs,  elemsets and  other parts  of  C<PETSc-FEM> are
passed via text-hashes, which are  then queried by the code via calles
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

Mario A. Storti E<lt>mstorti@intec.unl.edu.arE<gt>

=cut


