#!/usr/bin/perl

$case = 'vtube';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# GMV input file
if (1) {
    $rslt = "$case.state.good"; # petscfem result file
#    $rslt = "$case.ini.tmp"; # petscfem result file
    $rec = 0;
    $fields = 'nsc';
}

require "./data.pl";
require "$ENV{'HOME'}/PETSC/petscfem/tools/gmv.pl";
