#!/usr/bin/perl

system "octave -qH proc5.m";
$case = 'step3d';
$nod = "$case.defo_nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file
# $rslt = "outvector0.out.tmp"; # petscfem result file
# $rec = 'last';
# $fields = 'scalar';
# $fields = 'ns';

require '../../tools/gmv.pl';
