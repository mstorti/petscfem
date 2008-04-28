#!/usr/bin/perl

$case = 'stabi';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file

# $rslt = "stabi.state.tmp"; # petscfem result file
# $rslt = "stabi.ini.tmp"; # petscfem result file
$rslt = "outvector0.out.tmp"; # petscfem result file
$rec = 'last';

# $fields = 'scalar';
$fields = 'ns';

require '../../../tools/gmv.pl';
