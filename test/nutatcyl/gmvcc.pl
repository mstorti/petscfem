#!/usr/bin/perl

require 'math.pl';
# require './vaughn.pl';
$Omega = 3000/60*2*$PI;

$field_transf = sub {
    my ($u,$v,$w,$p) = @{$_[0]};
    my ($x,$y,$z) = @{$_[1]};
    my $uu = $u + $Omega * $y;
    my $vv = $v - $Omega * $x;
    return  $uu,$vv,$w,$p ;
};

$case = 'cylinder';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file

$nfic=2;			# eliminar ultimos $nfic nodos de la tabla de nodos
$nficcon=2;			# eliminar ultimos $nficcon nodos de la lista
				# de conectividades de cada elemento
$nficdof=2;			# eliminar ultimos $nficcon nodos de la lista
				# de conectividades de cada elemento


$rslt = "outvector.out.tmp";	# petscfem result file
$rec = 'last';
# $fields = 'scalar';
$fields = 'ns';

require "$ENV{'PETSCFEM_DIR'}/tools/gmv.pl";
