#!/usr/bin/perl

$case = 'vtube';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# GMV input file
if (1) {
    $rslt = "$case.state_0.tmp"; # petscfem result file
    $rec = 0;
    $fields = 'nsc';
}

require "./data.pl";
require "/u/mstorti/PETSC/petscfem.mario/tools/gmv.pl";
