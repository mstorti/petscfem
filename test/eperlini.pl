require "$ENV{'HOME'}/perl/math.pl";

sub heredoc {
    my $txt = shift();
    my $name = shift();
    $name = "tmp_file_$$.tmp" unless $name;
    open TMP,">$name";
    print TMP $txt;
    close TMP;
    print $name;
}

sub makeini {
    @U = @{shift()};
    $ndof = $#U;
    $nnod = shift();
    $name = shift();
    $noise = shift();
    print "name: $name\n";
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
# carl to eperl, or well it is the DEFAULT_VALUE. Pass the
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
    die "Couldn't find line \"$ident = <value>\"";
    
}

sub P {
    my $name = shift();
    my $val = shift();
    eval "\${\"$name\"} = $val";
    print "$name $val";
}

print <<'EOM';
# DON'T EDIT MANUALLY THIS FILE !!!!!!
# This files automatically generated by ePerl from 
# the corresponding `.epl' file. 
EOM
  
1;
