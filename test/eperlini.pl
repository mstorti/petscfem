#__INSERT_LICENSE__

use English;
require "$ENV{'PETSCFEM_DIR'}/tools/math.pl";
require "$ENV{'PETSCFEM_DIR'}/tools/utils.pl";
$NP = $ENV{'NP'};

# In `eperlini.pl' we have to check strings/number with `is_numeric()'
#   defined in 'perl/utils.pl' but this requires to update eperl in
#   spider. 

sub quote_string {
    my $value = shift();
    # This is how to know if a variable is a number or a string. 
    if ($value == 0 && $value !~ /^0$/) {
	return "\"$value\"";
    } else {
	return $value;
    }
}

## Gets a variable from the environment. If it was
## not set, set it to the default value
##
## usage: get_var_env(name,default)
sub get_var_env {
    my ($name,$def) = @_;
    ${$name} = $ENV{"$name"};
    ${$name} = $def unless ${$name};
}

sub get_var_env2 {
    my ($name,$def) = @_;
    ${$name} = (exists($ENV{"$name"}) ? $ENV{"$name"} : $def);
}

sub pr_auto {
    my $line = shift();
    my $begl = '';
    if ($line =~ /^(\W*)(\w)/) {
	$begl = $1;
	$line = $2.$POSTMATCH;
    }
    my @names = split " ",$line;
    my @text = ();
    foreach $name (@names) { push @text,"$name ",quote_string(${$name}),"\n"; }
    pop @text;
    print @text;
}

sub pr {
    my $line = shift();
    my @names = split " ",$line;
    my @text = ();
    foreach $name (@names) { push @text,"$name ",quote_string(${$name}),"\n"; }
    pop @text;
    print @text;
}

sub prc {
    my $line = shift();
    my @names = split " ",$line;
    my @text = ();
    foreach $name (@names) { push @text,"# ","$name ",quote_string(${$name}),"\n"; }
    pop @text;
    shift @text;
    print @text;
}

sub flag {
    $name = shift();
    print "$name ",(${$name} ? 1 : 0);
}

# usage: "doc_vals(VARS)", for instance "doc_vals(qw(Re,Ra,N,alpha))"
#
# includes legends of the form "# VAR VALUE" in the
# .depl file
sub doc_vals {
    my @vars = @_;
    print "#","-" x 20,"\n";
    foreach $var (@vars) {
	print "# \$$var = ${$var}\n" unless !defined($ {$var});
    }
    print "#","-" x 20,"\n";
}

sub heredoc {
    my $txt = shift();
    my $name = shift();
    $name = "tmp_file_$$.tmp" unless $name;
    open TMP,">$name";
    print TMP $txt;
    close TMP;
    print $name;
}

# usage: makeini(\@STATE,$nnod,$filename,$noise);
sub makeini {
    @U = @{shift()};
    $ndof = $#U;
    $nnod = shift();
    $name = shift();
    $noise = shift();
#    print "name: $name\n";
    if (!$name || !length($name)) {$name =  "tmp_file_$$.tmp";}
    open TMP,">$name";
    @UU = @U;
    for ($k=0; $k<$nnod; $k++) {
	if ($noise) {
	    @UU = @U;
  	    for ($kk=0; $kk<=$#UU; $kk++) {
  		$UU[$kk] *= (1+$noise*2*(rand()-1));
  	    }
	}
	print TMP "@UU\n";
    }
    close TMP;
    print $name;
}

# This is handy for setting options in the .epl data files. 
#
# usage: 
# <: option('OPTION_NAME',DEFAULT_VALUE) :>
# Then creates a line in the output of the form 
#
# OPTION_NAME VALUE
#
# where VALUE can be set through a `-d OPTION_NAME=VALUE' in the
# carll to eperl, or well it is the DEFAULT_VALUE. Pass the
# magic string __NULL__ if you want to pass a null string value. 
sub option {
    $option_name = shift();
    $option_default_val = shift();
    if ($ {$option_name}) {
	$option = $ {$option_name};
	$option="" if $option eq "__NULL__";
    } else {
	$option = $option_default_val;
    }
    print "$option_name $option";
}

#
# Same as before but only prints the value
#	
sub optval {
    $option_name = shift();
    $option_default_val = shift();
    if ($ {$option_name} ) {
	$option = $ {$option_name};
	$option="" if $option eq "__NULL__";
    } else {
	$option = $option_default_val;
    }
    return $option;
}

# usage readm(IDENT,FILE)
sub readm {
    my ($ident,$file)=@_;
    die "Couldn't open file $file\n" unless open MFILE,"$file";
    while (<MFILE>) {
	if (/^\s*$ident\s*=(\S*);/) {
	    close MFILE;
	    return $1;
	}
    }
    close MFILE;
    die "Couldn't find line \"$ident = <value>\"" unless $READM_RETURN_UNDEF;
    
}

sub import_vals {
    my ($file,$vars)=@_;
    die "Couldn't open file $file\n" unless open MFILE,"$file";
    my $read=1;
    print "# --- START READ FROM \"$file\" ----\n";
    while (<MFILE>) {
	if (/^.__ENDS_READING_VARS__$/) {
	    $read=0;
	    next;
	}
	if (/^.__START_READING_VARS__$/) {
	    $read=1;
	    next;
	}
	next if !$read;
	if (/^\s*(\w*)\s*=(\S*);/) {
	    $vars->{$1} = $2;
	    print "$1 $2\n"
	}
    }
    print "# --- END READ FROM \"$file\" ----\n";
}
    

sub readm_all {
    my @read_var=();
    my ($file)=@_;
    die "Couldn't open file $file\n" unless open MFILE,"$file";
    while (<MFILE>) {
	if (/^\s*(\w*)\s*=(\S*);/) {
	    push @read_var,$1;
	    $ {$1} = $2;
	}
    }
    close MFILE;
    return \@read_var;
}

sub readm_all_v {
    my @read_var=();
    my ($file)=@_;
    die "Couldn't open file $file\n" unless open MFILE,"$file";
    while (<MFILE>) {
	if (/^\s*(\w*)\s*=(\S*);/) {
	    push @read_var,$1;
	    $ {$1} = $2;
	}
    }
    close MFILE;
    return @read_var;
}

sub P {
    my $name = shift();
    my $val = shift();
    eval "\${\"$name\"} = $val";
    print "$name $val";
}

sub transcript {
    my $eperlfile = $ENV{'SCRIPT_SRC_PATH'};
#   If we wanted to transcript another file
#     ($octtmpfile,$eperlfile) = @_;
#     $eperlfile = $ENV{'SCRIPT_SRC_PATH'} unless $eperlfile;
    /`/;print <<EOM;
#
# Transcript of ePerl script: $eperlfile
# ===========================
#
EOM
/`/;
    open SCRIPT,$eperlfile;
    my $tr=0;
    while (<SCRIPT>) {
	$tr=0 if /^\#__END_TRANSCRIPT__$/;
	print "#  >$_" if $tr;
	$tr=1 if /^\#__TRANSCRIPT__$/;
    }
    close SCRIPT;

    return if $#_<0;
    my $octtmpfile = shift();
    if ($octtmpfile) {
	die "couldn't open $octtmpfile" unless open OCT,">$octtmpfile";
    }
    /`/; print <<EOM;
#
# [Eperlini library. "transcript" function.] Computed values:
# ===========================================================
EOM
/`/;
    foreach $v (@_) {
#	print "# \$$v: ${$v}\n";
	printf("# %10s: %s\n","\$$v",${$v});
	print OCT "$v = ${$v};\n" if $octtmpfile;
    }
    close OCT if $octtmpfile;
}

## Converts a string to an array of numbers (if posible),
## otherwise it leaves a string
## usage: $v = conv("1.2 2.3 3.5") -> returns [1.2,2.3,3.5]
##        $v = conv("A string") -> returns "A string"
sub conv {
    my $a0 = shift();
    chomp $a0;
    my $a = $a0; 
    my @values = ();
    while ($a =~ /\s*(\S*)/) {
	my $v = $1;
	my $rest = $POSTMATCH;
##	if ($v != 0 || $v =~ /^0(|.0*)(|e(|\+|\-)\d+)$/) {
	my $vv = getnum($v);
	if (defined $vv) { push @values,$vv; } 
	else { last; }
	$a = $rest;
    }
    if ($a =~ /^\s*$/) { return [@values]; }
    else { return $a; }
}

sub octave_export_vars {
    my ($file,@vars) = @_;
    my $exist =  -f $file;
    open O,"$file";
    print O <<'EOM' unless $exist && $file !~/^>>/;
# -*- mode: octave -*-  ## This is for Emacs
#  Automatically generated by 'octave_export_vars' (in eperlini.pl)
#  DO NOT EDIT MANUALLY THIS FILE !!!!
EOM
    for $name (@vars) {
	my $value = $ {$name};
	$value = conv($value);
	# This is how to know if a variable is a number or a string. 
	if (!ref($value)) {
	    print O "$name = \"$value\";\n"; # not number
	} elsif (@$value == 1) {
	    print O "$name = $value->[0];\n"; # number
	} else {
	    print O "$name = [",join(",",@$value),"];\n";
	}
    }
    close O;
}

# Added 'octave_string' that strips double quotes from
# arguments read from
# Octave data files. 
# usage octave_string(IDENTIFIERS...).
# e.g.: octave_string(qw(var1 var2 var3))
sub octave_string {
    my @args=@_;
    for $id (@args) {
	next unless defined(${$id});
	my $t=${$id};
	$t=~ s/^\s*\"//;
	$t=~ s/\"\s*$//;
	${$id} = $t;
    }
}

# usage: $array2 = aload("myfile")
#   returns an array ref to a 2D array 
#   after, gets elements with $array->[$i][$j]
sub aload {
    my $file = shift();
    die "can't open $file\n" unless open F,$file;
    my @data;
    while(<F>) { 
#	my @v = map { $_ = getnum($_); } split(" ",$_);
	next if /^\s*\#/;
	my @v = split(" ",$_);
	if (@v==1) {
	    push @data,$v[0]; 
	} else {
	    push @data,[@v]; 
	}
    }
    return \@data;
}

# usage: asave("myfile",$array2)
#   saves an array ref to a 2D array in file $file. 
#   Example: asave("myfile.dat",[[1,2,3],[4,5,6]]);
#   creates a file "myfile.dat" with rows "1 2 3" and "4 5 6". 
sub asave {
    my ($file,$array) = @_;
    die "can't open \"$file\"\n" unless open OUT,">$file";
    for my $r (@$array) {
	for my $e (@$r) {
	    print OUT $e;
	}
	print OUT "\n";
    }
    close OUT;
}

## Functions that help in writing hooks
## You simply call `add_hook(<hook-pair>);'
## When done, call `print_hook'
## For instance: 
##   <:add_hook('my_hook1 hook1'):>
##   ...
##   <:add_hook('my_hook2 hook2'):>
##   ...
##   <:add_hook('ns_dx_hook hook3'):>
##   ...
##   <:print_hook:>
sub add_hook {
    @hook_list = () unless defined @hook_list;
    push @hook_list,@_;
}

sub print_hook {
    my @t = ();
    if (@hook_list) {
	push @t,"hook_list ";
    }
    
    my $j=0;
    while (my $h = shift @hook_list) {	
	push @t,"     " if $j++;
	push @t,$h," \\\n"; 
    }
    pop @t;
    push @t,"\n";
    print @t;
}

print <<'EOM';
# DON'T EDIT MANUALLY THIS FILE !!!!!!
# This files automatically generated by ePerl from 
# the corresponding `.epl' file. 
EOM
  
1;
