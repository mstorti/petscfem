#!/usr/bin/perl

$case = 'cubcav';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file

## number of fictitious nodes
$nfic=0;			# eliminar ultimos $nfic nodos de la tabla de nodos
$nficcon=0;			# eliminar ultimos $nficcon nodos de la lista
				# de conectividades de cada elemento
$nficdof=0;			# eliminar ultimos $nficdof nodos de la lista
				# de conectividades de cada elemento
if (1) {
    # Load the last vector in $rslt
    $rslt = "outvector0.out.tmp"; # petscfem result file
    $rec='last';
    # $mat = "proc.tmp";
    # $mat = "ntype.tmp";
} else {
    # Load a specific vector
    $rec=20;
    $rslt = "state0.all.tmp";	# petscfem result file
}

require '../../tools/gmv.pl';
