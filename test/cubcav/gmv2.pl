#!/usr/bin/perl

$case = 'cubcav';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file
$rslt = "outvector0.out.tmp"; # petscfem result file
$rec = 'last';
$fields = 'scalar';

require '../../tools/gmv.pl';
