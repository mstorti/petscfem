# -*- perl -*-
#__INSERT_LICENSE__
#$Id: gmv.pl,v 1.8 2003/01/08 15:54:26 mstorti Exp $

if (! defined $fields) { $fields = 'ns'; }

open NOD,"$nod";
open GMV,">$gmv";

sub entropy_hook {
    my ($state) = @_;
    my $s = $state->[3]/$state->[0]**$ga;
    return $s;
}

sub enthalpy_hook {
    my ($state) = @_;
    my $h = $state->[3]/$state->[0]*$ga/($ga-1)+0.5*($state->[1]**2+$state->[2]**2);
    return $h;
}

sub print_rslt {
    my ($rslt,$nrslt,$nnod,$dim0,$dim1,$rec,$file) = @_;
    for (my $k=$dim0; $k<=$dim1; $k++) { 
	for (my $j=0; $j<$nnod; $j++) { 
	    print $file $rslt->[$j*$nrslt+$k],"\n";
	}
    }
}

@nodos = ();
$dim = 0;
while(<NOD>) {
    my @val = split " ";
    if ($#val>=0) {
	$dimm = $#val+1;
	$dim = $dimm if $dim==0;
	die "not consistent number of rows. Previously $dim ",
	"now $dimm.\n Recently read: \"$_\"" if $dimm!=$dim;
	push @nodos,@val;
    }
}
close NOD;

## remove fictitious nodes
for (my $k=0; $k<$nfic*$dim; $k++) { pop @nodos; }

$nnod = ($#nodos+1)/$dim;
$nnod_tot = ($#nodos+1)/$dim+$nfic;

/`/;
open GMV,">$gmv";
print GMV <<EOF;
gmvinput ascii
comments
  these are nodes
endcomm
nodes $nnod
EOF
/`/;

for (my $j=0; $j<$dim; $j++) {
    for (my $i=0; $i<$nnod; $i++) {
	print GMV $nodos[$dim*$i+$j],"\n";
    }
}

# GMV requires always three dimensions
for (my $j=$dim; $j<3; $j++) {
    for (my $i=0; $i<$nnod; $i++) {
	print GMV "0.\n";
    }
}

#undef $nodos;

open CON,"$con";
@cone=();
while(<CON>) {
    my @val = split " ";
    if ($#val>=0) {
	splice @val,$#val-$nficcon+1,$nficcon;
	$nelm = $#val+1;
	$nel = $nelm if $nel==0;
	die "not consistent number of rows. Previously $nel ",
	"now $nelm.\n Recently read: \"$_\"" if $nelm!=$nel;
	push @cone,@val;
    }
}
close CON;
$nelem = ($#cone+1)/$nel;
print "after reading connectivities...\n";

if ($dim==2 && $nel==4) {
    $shape = 'quad';
} elsif ($dim==3 && $nel==8) {
    $shape = 'hexa';
} elsif ($dim==3 && $nel==4) {
    $shape = 'tetra';
} else {
    die "unknown case dim $dim, nel $nel\n";
}

print GMV "cells $nelem\n";
for (my $ele=0; $ele<$nelem; $ele++) {
    print GMV "$shape $nel\n";
    for (my $i=0; $i<$nel; $i++) {
	print GMV ' ',$cone[$ele*$nel+$i];
    }
    print GMV "\n";
}

if ($rslt) {
    if ($rec eq 'last') {
	## count records in RSLT
	open RSLT,"$rslt";
	$nrec=0;
	while (<RSLT>) { $nrec++; }
	die "Not complete records in $rslt\n" unless $nrec % $nnod_tot ==0;
	$nrec = $nrec/$nnod_tot;
	$rec = $nrec;
    }

# Resultados
    open RSLT,"$rslt" || die "can't open result file $rslt\n";
    @rslt = ();			# result vector
    for (my $j=0; $j<$nnod_tot*($rec-1); $j++) { 
	print;
	<RSLT> || die "not enough records in file\n"; 
    }				# skip previous records

    for ($j=0; $j<$nnod; $j++) { 
	my $line = <RSLT>; 
	$line || die "end of archive detected\n";
	my @val = split " ",$line;
	splice @val,$#val-$nficdof+1,$nficcon;
	$nfield = scalar @val;
	if ($j==0) { $nrslt = $#val+1; }
	if (defined $field_transf) { 
	    my ($i1,$i2) = ($dim*$j,$dim*($j+1)-1);
	    my @xnod = @nodos[$i1..$i2];
	    @val = &$field_transf(\@val,\@xnod); 
	}
#	print "val: ",join(" ",@val),"\n";
	push @rslt,@val;
    }

    if ($fields eq 'ns') {
	print GMV "velocity 1\n";
	print_rslt(\@rslt,$nrslt,$nnod,0,$dim-1,0,GMV);
	for (my $j=$dim; $j<3; $j++) {
	    for (my $i=0; $i<$nnod; $i++) {
		print GMV "0.\n";
	    }
	}

	print GMV "variables\n";
	print GMV "pressure 1\n";
	print_rslt(\@rslt,$nrslt,$nnod,$dim,$dim,0,GMV);
	print GMV "endvars\n";

    } elsif ($fields eq 'nsc') {
	print GMV "velocity 1\n";
	print_rslt(\@rslt,$nrslt,$nnod,1,$dim,0,GMV);
	for (my $j=$dim; $j<3; $j++) {
	    for (my $i=0; $i<$nnod; $i++) {
		print GMV "0.\n";
	    }
	}

	print GMV "variables\n";
	print GMV "pressure 1\n";
	print_rslt(\@rslt,$nrslt,$nnod,$dim+1,$dim+1,0,GMV);
	print GMV "density 1\n";
	print_rslt(\@rslt,$nrslt,$nnod,0,0,0,GMV);
	if ($compute_s) {
	    my @s = ((0) x $nnod);
	    for (my $j=0; $j<$nnod; $j++) { 
		my $n0=$j*$nrslt;
		my $n1=($j+1)*$nrslt-1;
		my @state = @rslt[$n0..$n1];
		$s[$j] = entropy_hook(\@state);
#		print "entropy $s[$j]\n";
	    }
	    print GMV "entropy 1\n";
	    print_rslt(\@s,1,$nnod,0,0,0,GMV);
	}
	if ($compute_h) {
	    my @h = ((0) x $nnod);
	    for (my $j=0; $j<$nnod; $j++) { 
		my $n0=$j*$nrslt;
		my $n1=($j+1)*$nrslt-1;
		my @state = @rslt[$n0..$n1];
		$h[$j] = enthalpy_hook(\@state);
	    }
	    print GMV "enthalpy 1\n";
	    print_rslt(\@h,1,$nnod,0,0,0,GMV);
	}
	print GMV "endvars\n";

    } elsif ($fields eq 'scalar' || ! defined $fields) {
	
	print GMV "variables\n";
	for (my $f=1; $f<=$nfield; $f++) {
	    print GMV "phi$f $f\n";
	    print_rslt(\@rslt,$nrslt,$nnod,$f-1,$f-1,0,GMV);
	}
	print GMV "endvars\n";

    } else { die "unknown fields: \"$fields\"\n"; }

}

if ($mat) {
    @mat=();
    open MAT,"$mat";
    $nmat=0;
    while (<MAT>) { 
	if ($_ > $nmat) { $nmat = $_; }
	push @mat,$_+1; 
    }
    $nmat++;
    die "number of nodes do not match nnod $nnod, lines ",
    $#mat+1,"\n " if $#mat+1 != $nnod + $nfic;
    for ($k = 0; $k < $nfic; $k++) { pop @mat; }
    print GMV "material $nmat 1\n";
    for ($m=0; $m < $nmat-1; $m++) { print GMV "proc$m\n"; }
    print GMV "mixed\n";
    for $f (@mat) { print GMV "$f\n"; }
}

print GMV "endgmv\n";

close GMV;

print GMV "endgmv\n";

close GMV;
1;
