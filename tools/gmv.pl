# -*- perl -*-
#__INSERT_LICENSE__
#$Id: gmv.pl,v 1.4 2002/08/03 00:48:12 mstorti Exp $

if (! defined $fields) { $fields = 'ns'; }

open NOD,"$nod";
open GMV,">$gmv";

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
for (my $j=$dim; $j<3; $j++) {
    for (my $i=0; $i<$nnod; $i++) {
	print GMV "0.\n";
    }
}

undef $nodos;

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

if ($dim==2 && $nel==4) {
    $shape = 'quad';
} elsif ($dim==3 && $nel==8) {
    $shape = 'hexa';
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
    open RSLT,"$rslt";
    @rslt = ();			# result vector
    for (my $j=0; $j<$nnod_tot*($rec-1); $j++) { 
	<RSLT> || die "not enough records in file\n"; 
    }				# skip previous records

    for (my $j=0; $j<$nnod; $j++) { 
	my $line = <RSLT>; 
	$line || die "end of archive detected\n";
	my @val = split " ",$line;
	splice @val,$#val-$nficdof+1,$nficcon;
	$nfield = scalar @val;
	if ($j==0) { $nrslt = $#val+1; }
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
