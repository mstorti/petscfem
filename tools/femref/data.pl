require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl";

$N = 10; 
$nnod = ($N+1)**3; 
$nelem = $N**3*5;
$nfaces = $nelem*4;
$nfaces_ext = 2*6*$N**2;

@vars = qw(N);
octave_export_vars(">data.m.tmp",@vars);

system "octave -qH mkcube.m > ./mkcube.log";

$nnod_surf = count_lines("cube.surf-nod.tmp");

1;
