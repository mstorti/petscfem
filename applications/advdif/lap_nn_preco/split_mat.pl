#!/usr/bin/perl
#$Id: split_mat.pl,v 1.1.2.1 2004/03/02 17:21:59 mstorti Exp $
$nod = "nodos.dat";		# nodes
$seg = "segmentos.dat";		# segmentos
$tri = "triangulos.dat";		# triangs
$nodos2int = "n2int.dat";        # nodos para interpolar
$aq_con = "aquifer.con.tmp";        # aquif conect
$str_con = "stream.con.tmp";        # stream conect
$str_loss_con = "stream_loss.con.tmp";        # stream_loss conect
$fixa = "fase_I.fixa.tmp"; #fijaciones

open NOD,"$nod";
#pattern= {nnod, tipo, #, x, y , h, f} cambiar esto si cambia el archivo de entrada
#necesito armar el mapeo entre los indices de los nodos que despues seran nodos del stream!!!!
$count_stn = 0;
$count_aqbn = 0;
@nodos = ();
@rep = ();
$dim = 2;
$line = <NOD>;
@fields = split " ",$line;
$nnod = shift @fields;
$naq = $nnod; # #aquifer nodes
%a2s = ();
# pp == 1 : salado (frontera fase I)
# pp == 3 : naciente salado
# pp == 4 : rios internos
# pp == 5 : desembocadura rios internos a la frontera salado
# pp == 12: nacientes rios internos
@nac = ();
$rec_chan = 5;
$rec_sal = 5;
$u_chan = 0;
$u_sal = 0;

$dof_u = 1;
$dof_h = 2;
$c_nac = 0;
%aq_b_n = ();
for (my $k=0; $k<$nnod; $k++) { 
    $line = <NOD>;
    my @ligne = split " ",$line;
    if (@ligne>=0) {
	push @nodos, @ligne[1..$dim];
	my $pp = $ligne[4];
	if ($pp == 1 || $pp == 3 || $pp == 4 || $pp == 5 || $pp == 12) {
	    $count_stn++;
	    push @rep, @ligne[1..$dim];
            # aca armo el map
	    $a2s{$ligne[0]} = $naq+$count_stn;
	}
	if ($pp == 2){
	    $count_aqbn++;
	    $aq_b_n{$ligne[0]} = $count_aqbn;
	}
    }
}
close NOD;
#$count_aqbn = @aq_b_n;

open NOD, "$nod";
for (my $k=0; $k<$nnod; $k++) { 
    $line = <NOD>;
    my @ligne = split " ",$line;
    if (@ligne>=0) {
	my $pp = $ligne[4];
	if ($pp == 3 && $ligne[0] != 0){
	    $c_nac++;
	    push @nac, $a2s{$ligne[0]},$dof_u,$u_sal;
	}
	if ($pp == 12) {
	    $c_nac++;
	    push @nac, $a2s{$ligne[0]},$dof_u,$u_chan;
	}
    }
}
close NOD;

$nnrep = @rep;
$nstr = $nnrep/2; # #river nodes

if (0){
open PP,">pp.nod";
for(my $k=0;$k<$nstr;$k++){
    print PP $rep[2*$k], ' ' ,$rep[2*$k+1],"\n";
}
close PP;
}
#pongo los rep(==rios) al final de nodos.. + el nodo dummy para el mult lag 
for(my $k=0;$k<$nnrep;$k++){
    push @nodos, $rep[$k];
}
undef $nnod;
$nnod = (@nodos/$dim);
open N2I,">$nodos2int";
for (my $j=0; $j<$nnod; $j++) {
    for (my $i=0; $i<$dim; $i++) {
	print N2I $nodos[$dim*$j+$i],"\t\t";
    }
    print N2I "\n";
}
close N2I;

open CON,"$tri";
$vertex = 3; 
@cone=();
$line = <CON>;
@fields = split " ",$line;
$ncone = shift @fields;
for (my $k=0; $k<$ncone; $k++) { 
    $line = <CON>;
    my @ligne = split " ",$line;
    if (@ligne>=0) {
	push @cone, @ligne[1..$vertex];
    }
}
$nelm = @cone/$vertex;
close CON;

$ele0 = "aquifer.con0.tmp";
open EL0, ">$ele0";
open ELE,">$aq_con";
for (my $ele=0; $ele<$nelm; $ele++) {
    for (my $i=0; $i<$vertex; $i++) {
	print ELE ' ',$cone[$ele*$vertex+$i]+1; #se suma 1 para llevar a base 1
	print EL0 ' ',$cone[$ele*$vertex+$i]; 
    }
    print ELE "\n";
    print EL0 "\n";
}
close EL0;
close ELE;


open SEG,"$seg";
@fields = ();
$vertex = 2; 
@cone = ();
$line = <SEG>;
@fields = split " ",$line;
$ncone = shift @fields;
for (my $k=0; $k<$ncone; $k++) { 
    $line = <SEG>;
    my @ligne = split " ",$line;
    if (@ligne >= 0 && (!exists $aq_b_n{$ligne[1]} && !exists $aq_b_n{$ligne[$vertex]})) {
	push @cone, @ligne[1..$vertex];
    }
}
$nelm = @cone/$vertex;
close SEG;

open STR,">$str_con";
open STL,">$str_loss_con";
for (my $ele=0; $ele<$nelm; $ele++) {
    my $pp=$cone[$ele*$vertex];
    my $pp1=$cone[$ele*$vertex+1];
    if (exists $a2s{$pp} && exists $a2s{$pp1}) {
	print STR ' ',$a2s{$pp},' ',$a2s{$pp1};
	print STL ' ',$a2s{$pp},' ',$a2s{$pp1}, ' ',$pp+1,' ',$pp1+1;
	#no se suma 1 en los rios para llevar a base 1 ya que se hizo en el hash ($naq) 
    }
    print STR "\n";
    print STL "\n";
}
close STR;
close STL;
open FIX, ">$fixa";
$n_fix = 3;
for (my $i=0; $i<$c_nac;$i++){
    for (my $j=0; $j<$n_fix; $j++) {
	print FIX ' ',$nac[$n_fix*$i+$j];
    }
    print FIX "\n";
}
close FIX;


if (1){
    $n2i = "n2int.dat";
#    $z_int = "z_int_n2int.nod";		# nodes
    $z_int = "fase_I.nod_fict_int.tmp";		# nodes
    $femnode = "fase_I.nod_fict.tmp";		# segmentos
    $ini = "fase_I.ini.tmp";
    open N2I,"$n2i";
    open ZIN,"$z_int";
    open FEM,">$femnode";
    open INI, ">$ini";
    $eta_aq = 45.0;
    $h_b = 10.0;
    $aq_ini = 30.0;
    $str_ini = 5.0;
    $cc = 0;
    while (<ZIN>) {
	$cc++;
	my $z = $_;
	$line = <N2I>;
	my @ligne = split " ",$line;
	if (@ligne >= 0) { 
	    if ($cc < $naq+1) {
		print FEM $ligne[0],' ',$ligne[1], ' ',$z-$eta_aq, "\n";
		print INI 0.0 , ' ', $z-$eta_aq+$aq_ini, "\n"; 
	    } else {
		print FEM $ligne[0],' ',$ligne[1], ' ',$z-$h_b, "\n";
		print INI $u_chan, ' ', $str_ini, "\n";
	    }
	}
    }
    $x = 5428710;
    $y = 6496698;
    $eta_dummy = 0;
    print FEM $x, ' ', $y,' ', $eta_dummy, "\n";
    print INI 0, ' ', 0, "\n";
    close N2I;
    close ZIN;
    close FEM;
    close INI;
}

$rest = "salado.rest.tmp";
open RES, ">$rest";
$nndummy = $naq+$nstr+1;
print RES $naq+1, ' ', $naq+2, ' ', $naq+3, ' ', $naq+4, ' ', $nndummy, ' ', 2.000001 , ' ', 5.0000001, ' ', 1.0000000001, "\n"; 
close RES;
