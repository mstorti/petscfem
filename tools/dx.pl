# -*- perl -*-
#__INSERT_LICENSE__
#$Id: dx.pl,v 1.3 2003/01/29 18:06:27 mstorti Exp $

sub connect_0 {
    my ($file,$file_out,$hook) = @_;
    if (!defined $file_out) {
	if (-f "$file~") { unlink "$file~"; }
	rename $file,"$file~";
	$file_out = $file;
	$file = "$file~";
    }
    open FILE,$file;
    open FILEO,">$file_out";
    my $lines=0;
    while (<FILE>) { 
	my @l = split " ";
	&$hook(\@l) if defined $hook;
	for $f (@l) { print FILEO $f-1," "; }
	print FILEO "\n";
	$lines++;
    }
    close FILE;
    close FILEO;
}

sub quad_hook {
    my $r = shift();
    my $n3 = $r->[3];
    $r->[3] = $r->[2];
    $r->[2] = $n3;
}

sub ns_extract {
    my ($fin,$uf,$pf) = @_;
    open FIN,"$fin";
    open FU,">$uf";
    open FP,">$pf";
    my $ndim;
    while (<FIN>) {
	my @l = split " ";
	if (!defined $ndim) { $ndim = @l-1; }
	for ($j=0; $j<$ndim; $j++) { print FU "$l[$j] "; }
	print FU "\n";
	print FP "$l[$ndim]\n";
    }
    close FIN;
    close FU;
    close FP;
}

sub nsc_extract {
    my ($fin,$uf,$pf,$rhof) = @_;
    open FIN,"$fin";
    open FU,">$uf";
    open FP,">$pf";
    open FRHO,">$rhof";
    my $ndim;
    while (<FIN>) {
	my @l = split " ";
	if (!defined $ndim) { $ndim = @l-2; }
	print FRHO "$l[0]\n";
	for ($j=0; $j<$ndim; $j++) { print FU "$l[1+$j] "; }
	print FU "\n";
	print FP "$l[$ndim+1]\n";
    }
    close FIN;
    close FU;
    close FP;
    close FRHO;
}

sub nsther_extract {
    my ($fin,$uf,$pf,$Tf) = @_;
    open FIN,"$fin";
    open FU,">$uf";
    open FP,">$pf";
    open FT,">$Tf";
    my $ndim;
    while (<FIN>) {
	my @l = split " ";
	if (!defined $ndim) { $ndim = @l-2; }
	for ($j=0; $j<$ndim; $j++) { print FU "$l[$j] "; }
	print FP "$l[$ndim]\n";
	print FT "$l[$ndim+1]\n";
    }
    close FIN;
    close FU;
    close FP;
    close FT;
}

sub count_lines {
    my $file = shift();
    open FILE,$file;
    my $lines=0;
    while (<FILE>) { $lines++; }
    return $lines;
}

sub prism_split_dx {
    my ($prism_c,$tetra_c) = @_;
    open PRISM,"$prism_c";
    open TETRA,">$tetra_c";
    my $prisms=0;
    while (<PRISM>) {
	next if /^\#/;
	$prisms++;
	my @nodes = split " ";
	for (my $k=0; $k<6; $k++) { $nodes[$k] -= 1; }
	die "numero de nodos debe ser 6!!" unless @nodes==6;
	print TETRA "$nodes[0] $nodes[1] $nodes[2] $nodes[3]\n";
	print TETRA "$nodes[4] $nodes[3] $nodes[5] $nodes[1]\n";
	print TETRA "$nodes[1] $nodes[5] $nodes[2] $nodes[3]\n";
    }
    close TETRA;
    close PRISM;
}

1;
