#!/usr/bin/perl

require "utils.pl";
$case = 'cylin';
$nod = "$case.nod.tmp";	# nodes
$con = "$case.con.tmp";	# connectivities
$gmv = "$case.gmv.tmp";	# MV input file
# $rslt = "cylin.ini.tmp"; # petscfem result file
$rslt = "cylin.state.tmp"; # petscfem result file
$rec = 'last';
# eliminar ultimos $nfic nodos de la tabla de nodos
$nfic = count_lines('ext.coupling_nodes.tmp');
# $fields = 'scalar';
$fields = 'ns';

require '../../tools/gmv.pl';
