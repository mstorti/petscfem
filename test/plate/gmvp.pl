#!/usr/bin/perl

$case = 'ext';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file
# $rslt = "cylin.ini.tmp"; # petscfem result file
$rslt = "ext.state.tmp"; # petscfem result file
$rec = 'last';
$fields = 'scalar';
# $fields = 'ns';

require '../../tools/gmv.pl';
