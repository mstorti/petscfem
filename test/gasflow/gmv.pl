#!/usr/bin/perl

$case = 'nozzle';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# GMV input file
$rslt = "$case.state.tmp"; # petscfem result file
$rec = 0;
$fields = 'nsc';

require "./data.pl";
require "$ENV{'HOME'}/PETSC/petscfem/tools/gmv.pl";
