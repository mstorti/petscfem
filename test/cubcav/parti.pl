$nsubdo = 30;
for ($s = 0; $s<$nsubdo; $s++) {
    $nelem = count_lines("icone$s.dat");
    print <<EOF;
#================================================================
object "s$s" class array type int rank 1 shape 4
	items $nelem data file "icone$s.dat"
attribute "element type" string "tetrahedra"
attribute "ref" string "positions"

object "colors$s" class array type float rank 1 shape 3 items $nelem
           data file "colors$s.dat"
attribute "dep" string "connections"

object "subdo$s" class field 
component "positions" value "xnod"
component "connections" value "s$s"
component "colors" value "colors$s"

EOF
}

$ngroup = 6;
$sg = ceil($nsubdo/$ngroup);
for ($g = 0; $g<$ngroup; $g++) {
    $s1 = $g*$sg;
    $s2 = $s1+$sg;
    $s2 = $nsubdo if $s2>$nsubdo;
    print "object \"sg$g\" class group\n";
    for ($j=$s1; $j<$s2; $j++) {
	print "member \"subdo$j\" value \"subdo$j\"\n";
    }
    print "\n";
}

print "object \"group\" class group\n";
for ($g = 0; $g<$ngroup; $g++) {
    print "member \"sg$g\" value \"sg$g\"\n";
}
print "\n";
print "end\n";

1;
