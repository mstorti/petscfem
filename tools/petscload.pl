#!/usr/bin/perl
#  [-begin-license-]
#  This  file  belongs  to  the  PETSc  -  FEM  package,  a  library  and
#  application  suite oriented  to  the Finite  Element  Method based  on
#  PETSc.   Copyright  (C)  1999-2001,  Mario  Alberto  Storti,  Norberto
#  Marcelo Nigro, Centro Internacional de Metodos Numericos en Ingenieria
#  (CIMEC-Argentina),  Universidad Nacional del  Litoral (UNL-Argentina),
#  Consejo   Nacional   de   Investigaciones   Cientificas   y   Tecnicas
#  (UNL-Argentina).
#  
#  This program is  free software; you can redistribute  it and/or modify
#  it under the  terms of the GNU General Public  License as published by
#  the Free Software Foundation; either  version 2 of the License, or (at
#  your option) any later version.
#  
#  This program  is distributed in the  hope that it will  be useful, but
#  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
#  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
#  General Public License for more details.
#  
#  You  should have received  a copy  of the  GNU General  Public License
#  along  with  this  program;  if   not,  write  to  the  Free  Software
#  Foundation, Inc.,  59 Temple Place, Suite 330,  Boston, MA 02111-1307,
#  USA.
#  [-end-license-]



$in_data=0;
$name="";
$basename="tmp_petsc_data";
@namelist=();

while(<>) {
    if (/^(\S*) = \[/) {
	$in_data=1;
	$name=$1;
	open TMP,">$basename._";
    } elsif (/\]/) {
	$in_data=0;
	close TMP;
	rename "$basename._","$basename.$name";
	push @namelist,$name;
    } elsif ($in_data) {
	print TMP;
    } elsif (/\s*(\S*)\s*=\s*spconvert/) {
	$nameo=$name;
	$name=$1;
	pop @namelist;
	push @namelist,$name;
	rename "$basename.$nameo","$basename.$name";
    } elsif (/^%/) {
    } elsif (/^zzz = zeros/) {
    } else {
	print "????? -> <$_>\n";
    }
}
	
open TMP,">${basename}_script.m";
foreach $name (@namelist) {
    print TMP "load $basename.$name; $name=$basename; ",
    "unlink(\"$basename.$name\"); clear $basename;\n";
}
close TMP;
