
require '../../tools/dx.pl';

$xnod = "vtube.nod-2d.tmp";
$icone = "vtube.con0-2d.tmp";
$nnod = count_lines($xnod);
$nelem = count_lines($icone);
nsc_extract("./vtube.state.tmp","u.dat.tmp","p.dat.tmp","rho.dat.tmp");

1;
