$ok=1;
while (<STDIN>) {
    if (! /   filen(\d)*.dat:(\d*)   <line(\d*) in filen(\d*)>/ 
	|| $1!=$4 || $2!=$3) {
	print "bad line : \"$_\"\n";
	$ok=0;
    }
#    print " $1 $2 $3 $4\n";
}
print "file ok? $ok\n";

