
require '../../tools/dx.pl';

$dir = "./" unless defined $dir;
$dir .= "/" unless $dir =~ m|/$|;
$state = "./vtube.state.tmp";
if (defined $step && $step>=0) {
    if (-f "${dir}vtube.state_$step.tmp.gz") {
	print "# unpacking \"${dir}vtube.state_$step.tmp.gz\"\n";
	system "gunzip -c ${dir}vtube.state_$step.tmp.gz > ./vtube.tempo.tmp";
	$state = "./vtube.tempo.tmp";
    } elsif (-f "${dir}vtube.state_$step.tmp") {
	$state = "${dir}vtube.state_$step.tmp";
    } else { die "can't find \"${dir}vtube.state_$step.tmp(|.gz)\n"; }
}
print "# using state: $state\n";
$xnod = "vtube.nod-2d.tmp";
$icone = "vtube.con0-2d.tmp";
$nnod = count_lines($xnod);
$nelem = count_lines($icone);
nsc_extract($state,"u.dat.tmp","p.dat.tmp","rho.dat.tmp");

1;
