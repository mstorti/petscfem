#__INSERT_LICENSE__

require "$ENV{'PETSCFEM_DIR'}/tools/math.pl";
$NP = $ENV{'NP'};

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
## not set, set it to de default value
##
## usage: get_var_env(name,default)
sub get_var_env {
    my ($name,$def) = @_;
    ${$name} = $ENV{"$name"};
    ${$name} = $def unless ${$name};
}

sub pr {
    $name = shift();
    print "$name ",quote_string(${$name});
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
	print "# \$$v: ${$v}\n";
	print OCT "$v = ${$v};\n" if $octtmpfile;
    }
    close OCT if $octtmpfile;
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
	# This is how to know if a variable is a number or a string. 
	if ($value == 0 && $value !~ /^0$/) {
	    print O "$name = \"$value\";\n"; # string
	} else {
	    print O "$name = $value;\n"; # number
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

print <<'EOM';
# DON'T EDIT MANUALLY THIS FILE !!!!!!
# This files automatically generated by ePerl from 
# the corresponding `.epl' file. 
EOM
  
1;
