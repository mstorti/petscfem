#!/usr/bin/perl
#__INSERT_LICENSE__

require '../tools/myexpect.pl';

use Getopt::Std;
getopts("hond");

if ($opt_h) {
    print <<'EOM';
usage: ` $ runtests.pl [OPTIONS]'

OPTIONS: -h print help
         -o complain when can't open files
         -n report only 'Not OK' tests
         -d debug
EOM
    exit;
}

$COMPLAIN_ON_CANT_OPEN = $opt_o;

begin_section('All tests');

begin_section('FastMat2');

##------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
/'/; # to disable auto-filling an indenting in Emacs
expect("fastmat2/fastmat2a.out.tmp","FastMat2 library",
        read_file("fastmat2/fastmat2a.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fastmat2/fastmat2.out.tmp","FastMat2 library with caches",
                         read_file("fastmat2/fastmat2.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fastmat2/fastmat2c.out.tmp","FastMat2 library. FEM problem.",<<'EOT');
total_area
Total area OK\? : YES
Summary of operation counts:
     get:  24
     put:  12
     sum:  4
     mult: 0
     div:  0
     abs:  0
     fun:  0
CPU:
rate:
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Laplace boundary conditions');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/
$sector_test = <<'EOT';
-6\.787.*e-02  6\.787.*e-02  
-1\.359.*e-01  1\.359.*e-01  
-2\.043.*e-01  2\.043.*e-01  
-2\.732.*e-01  2\.732.*e-01  
-3\.427.*e-01  3\.427.*e-01  
-4\.132.*e-01  4\.132.*e-01  
-4\.846.*e-01  4\.846.*e-01  
-5\.573.*e-01  5\.573.*e-01  
-6\.314.*e-01  6\.314.*e-01  
-7\.071.*e-01  7\.071.*e-01  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sector/sector0.sal","Periodic boundary conditions (1)",$sector_test);

expect("lap_per.sal","Periodic boundary conditions (2)",<<'EOT');
-8\.3417.*-02  8\.3417.*-02  
-1\.6401.*-01  1\.6401.*-01  
-2\.3928.*-01  2\.3928.*-01  
-3\.0721.*-01  3\.0721.*-01  
-3\.6642.*-01  3\.6642.*-01  
-4\.1628.*-01  4\.1628.*-01  
-4\.5728.*-01  4\.5728.*-01  
-4\.9418.*-01  4\.9418.*-01  
-5\.6250.*-01  5\.6250.*-01  
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Oscplate tests. Time dep. b.c.s and N.S.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Here we check only the $y$ components of velocity, since the others are
# rather spurious. 
expect("oscplate1.sal","Time dependent boundary conditions",<<'EOT');
4.76.*-01  
1.53.*-01  
2.38.*-03
-4.38.*-02 
-4.25.*-02 
-2.74.*-02 
-1.36.*-02 
-5.14.*-03 
-1.33.*-03 
EOT

#  5.100.*e-01
#  1.819.*e-01
#  1.980.*e-02
#  -3.595.*e-02
#  -4.028.*e-02
#  -2.769.*e-02
#  -1.451.*e-02
#  -5.920.*e-03
#  -1.748.*e-03

#   5\.1008.*e-01 
#   1\.8209.*e-01 
#   1\.9983.*e-02 
#  -3\.5744.*e-02 
#  -4\.0064.*e-02 
#  -2\.7483.*e-02
#  -1\.4343.*e-02
#  -5\.7924.*e-03
#  -1\.6814.*e-03


#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Same as before, but here the movemente of the plate is of
# saw-teeth wave-form. ($u=\pm constant$). 
expect("oscplate2o.sal","Time dependent boundary conditions (2)",<<'EOT');
size is OK\? 1
max\(u\)<tol OK\? 1, 
max\(p\)<tol OK\? 1, 

v en la pared, primer periodo
  -0.50000
  -0.50000
  -0.50000
  -0.50000
  -0.50000
  -0.50000
  -0.50000
  -0.50000
   0.50000
   0.50000
   0.50000
   0.50000
   0.50000
   0.50000
   0.50000
  -0.50000

v en x=0.5, primer periodo

   -7.2990e-04
   -3.4210e-03
   -8.7188e-03
   -1.6379e-02
   -2.5714e-02
   -3.6017e-02
   -4.6734e-02
   -5.7488e-02
   -6.6576e-02
   -7.1393e-02
   -7.0570e-02
   -6.4555e-02
   -5.4713e-02
   -4.2461e-02
   -2.8916e-02
   -1.6314e-02

v en la pared
  -0.50000
   0.50000
  -0.50000
  -0.50000
   0.50000
   0.00000
   0.00000
   0.00000

v en x=0.5
    2.3506e-02
   -2.1951e-02
    4.6242e-02
    1.1459e-02
    2.5354e-03
    5.6099e-04
    1.2413e-04
    2.7464e-05
    6.0768e-06
    1.3452e-06
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Oscillating plate with sinusoidal oscillation
# Verifies quadratic convergence of the
# Crank Nicholson method. Runs a case with N=16,
# 32 and 128 time steps. The solution for N=128 is
# considered the exact, and it is checked in Octave that
# norm(u32-u128)/norm(u16-u128) < 0.30 (it should be <0.25). 
#
expect("oscplate3o.sal","Time dep. b.c./Crank-Nicholson ",<<'EOT');
norm\(u32-u128\)/norm\(u16-u128\) < 0.26 OK \? > 1 
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Oscillating plate with sinusoidal oscillation
# periodic function interpolated with splines - "sin" case
expect("oscplate4o.sal","Time dep. b.c., spline_periodic function/sin",<<'EOT');
error = .*,  < .* OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Oscillating plate with sinusoidal oscillation
# periodic function interpolated with splines - "cos" case
expect("oscplate4co.sal","Time dep. b.c., spline_periodic function/cos",<<'EOT');
error = .*,  < .* OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Oscillating plate with sinusoidal oscillation
# periodic function interpolated with splines - "sin shift 30 dg" case
expect("oscplate4_30dego.sal",
"Time dep. b.c., spline_periodic function/sin shift",<<'EOT');
error = .*,  < .* OK \? 1
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Plano tests on advdif');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_local.sal","Adv/Sh.Water FastMat2/local_time_step",<<'EOT');
consistent_supg_matrix -> 0
__REWIND__
Courant -> 0.6
__REWIND__
local_time_step -> 1
__REWIND__
type *"volume_shallowfm2"
time_step 1, time: 0.6, res = 2.5634.*e-02
time_step 2, time: 1.2, res = 2.5538.*e-02
time_step 3, time: 1.8, res = 2.5456.*e-02
time_step 4, time: 2.4, res = 2.5386.*e-02
time_step 5, time: 3, res = 2.5327.*e-02
EOT

#old_results
<<PEPE;
time_step 1, time: 0.6, res = 2.5634
time_step 2, time: 1.2, res = 2.5764
time_step 3, time: 1.8, res = 2.5908
time_step 4, time: 2.4, res = 2.6066
time_step 5, time: 3, res = 2.6234
PEPE

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_auto.sal","Adv/Sh.Water FastMat2/auto_time_step",<<'EOT');
consistent_supg_matrix -> 0
__REWIND__
Courant -> 0.6
__REWIND__
auto_time_step -> 1
__REWIND__
local_time_step -> 0
type *"volume_shallowfm2"
time_step 1, time: 0.0196769, res = 2.56344.*e-02
time_step 2, time: 0.0393532, res = 2.55382.*e-02
time_step 3, time: 0.0590255, res = 2.54556.*e-02
time_step 4, time: 0.0786908, res = 2.53853.*e-02
time_step 5, time: 0.0983473, res = 2.53256.*e-02
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_local_nm.sal",
       "Adv/Sh.Water Newmat/local_time_step",<<'EOT');
consistent_supg_matrix -> 0
__REWIND__
Courant -> 0.6
__REWIND__
local_time_step -> 1
__REWIND__
type *"volume_shallow"
time_step 1, time: 0.6, res = 2.5634.*e-02
time_step 2, time: 1.2, res = 2.5538.*e-02
time_step 3, time: 1.8, res = 2.5456.*e-02
time_step 4, time: 2.4, res = 2.5386.*e-02
time_step 5, time: 3, res = 2.5327.*e-02
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_cons.sal",
       "Adv/Sh.Water FastMat2/cons_supg",<<'EOT');
consistent_supg_matrix -> 1
__REWIND__
auto_time_step -> 0
__REWIND__
local_time_step -> 0
__REWIND__
type *"volume_shallowfm2"
time_step 1, time: 0.1, res = 2.56344.*-02
time_step 2, time: 0.2, res = 2.96917.*-02
time_step 3, time: 0.3, res = 3.52613.*-02
time_step 4, time: 0.4, res = 4.27404.*-02
time_step 5, time: 0.5, res = 5.30517.*-02
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_cons_auto.sal",
       "Adv/Sh.Water FastMat2/cons_supg/auto",<<'EOT');
 -- Global_options: 
auto_time_step -> 1
__REWIND__
type *"volume_shallowfm2"
time_step 1, time: 0.00327948, res = 2.563449.*-02
time_step 2, time: 0.00655894, res = 2.563137.*-02
time_step 3, time: 0.00983839, res = 2.56286.*-02
time_step 4, time: 0.0131178, res = 2.562627.*-02
time_step 5, time: 0.0163972, res = 2.56242.*-02
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_local_weak.sal",
       "Adv/Sh.Water FastMat2/local_step/auto",<<'EOT');
 -- Global_options: 
__REWIND__
auto_time_step -> 0
__REWIND__
type *"volume_shallow"
time_step 1, time: 0.6, res = .*e-(15|16|17)
time_step 2, time: 1.2, res = .*e-(15|16|17)
time_step 3, time: 1.8, res = .*e-(15|16|17)
time_step 4, time: 2.4, res = .*e-(15|16|17)
time_step 5, time: 3, res = .*e-(15|16|17)
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("plano_fm2_weak.sal",
       "Adv/Sh.Water FastMat2/weak/bcconv_adv_fm2",<<'EOT');
 -- Global_options: 
__REWIND__
auto_time_step -> 0
__REWIND__
type *"volume_shallowfm2"
time_step 1, time: 0.6, res = .*e-(15|16|17)
time_step 2, time: 1.2, res = .*e-(15|16|17)
time_step 3, time: 1.8, res = .*e-(15|16|17)
time_step 4, time: 2.4, res = .*e-(15|16|17)
time_step 5, time: 3, res = .*e-(15|16|17)
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Navier-Stokes tests.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.re1000.verif.tmp",
       "Square cavity, Re1000, N=20",<<'EOT');
Weak form 0. error < tol OK \? 1
Weak form 1. error < tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.fs.verif.tmp",
       "Square cavity with fract. step, Re400, N=20",<<'EOT');
Square cavity at Re=400. Error < tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.zwproc.tmp",
       "Sqcav 0weight processor",<<'EOT');
Sq. Cavity with 0weight proc. error < tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.zwproc2.tmp",
       "Sqcav 0weight processor (iisd_subpart=2)",<<'EOT');
__EXACT_MATCH__
Sq. Cavity with 0weight proc (sbprt=2). 
error < tol OK ? 1 
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nutatcyl/checknutcyl.verif.tmp",
       "Nutating cylinder",<<'EOT');
Roll moment OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nutatcyl/rigid.verif.tmp",
       "Viscous force intgrator on rigid cylinder",<<'EOT');
Mz\(rigid\) OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nutatcyl/strip.verif.tmp",
       "Viscous force intgrator on linear 3D strip",<<'EOT');
Fx \[constant acceleration\] OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("mmove/step.verif.tmp",
       "Max volume ratio large deformation mesh mov.",<<'EOT');
Mesh OK \(all areas >0\) \? 1
Vol ratio < max allowed OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("mmove/step3d.verif.tmp",
       "Max volume ratio large deformation mesh mov.",<<'EOT');
All tetra volumes > 0 OK \? 1
Max/min ratio < max_ratio OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Generic load for NS elemset');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nsgenload/verif.scalar.tmp",
       "Generic load film elemset (scalar) ",<<'EOT');
Error < tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nsgenload/verif.diag.tmp",
       "Generic load film elemset (diag) ",<<'EOT');
Error < tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nsgenload/verif.full.tmp",
       "Generic load film elemset (full) ",<<'EOT');
Error < tol OK \? 1
EOT

end_section();

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Misc tests.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.iisd.verif.np1.tmp",
       "Square cavity, IISD part. in 1 proc.",read_file("sqcav/sqcav.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.g_body.verif.tmp",
       "Square cavity, test G_body term.",read_file("sqcav/gbody.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/check.iisd.verif.np2.tmp",
       "Square cavity, IISD part. in 2 proc.",read_file("sqcav/sqcav.ans.np2.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("graph/output.graph.tmp","Graph partitioning",
               read_file("graph/answer.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fastmat2/tfastvec.out.tmp",
       "FastVector library",<<'EOT');
resized  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  
resized  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  
with doubles  3.3  3.3  3.3  3.3  3.3  3.3  3.3  3.3  3.3  3.3  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("tidmap.sal",
       "Idmap class",<<'EOT');
testing for size of idmap n=5, number of operations 100
before permutations
checking consistency by rows: 
 ------------> OK\? : YES
checking consistency by cols: 
 ------------> OK\? : YES
Row permutations: maximum error: 
OK\? : YES
After 
checking consistency by rows: 
 ------------> OK\? : YES
checking consistency by cols: 
 ------------> OK\? : YES
Column permutations: maximum error: 
OK\? : YES
After set_elem\(\): maximum error:
OK\? : YES
After get_val\(\): maximum error: 
OK\? : YES
matrix for solve\(\)
After permuting and setting random: maximum error: 
OK\? : YES
After matrix solving: maximum error: 
OK\? : YES
testing for size of idmap n=.*, number of operations .*
checking consistency by rows: 
 ------------> OK\? : YES
checking consistency by cols: 
 ------------> OK\? : YES
Row permutations: maximum error: 
OK\? : YES
checking consistency by rows: 
 ------------> OK\? : YES
checking consistency by cols: 
 ------------> OK\? : YES
Column permutations: maximum error: 
OK\? : YES
checking consistency by rows: 
 ------------> OK\? : YES
checking consistency by cols: 
 ------------> OK\? : YES
Row rotations: maximum error: 
OK\? : YES
checking consistency by rows: 
 ------------> OK\? : YES
checking consistency by cols: 
 ------------> OK\? : YES
Column rotations: maximum error: 
OK\? : YES
After set_elem\(\): maximum error: 
OK\? : YES
After get_val\(\): maximum error: 
OK\? : YES
After permuting and setting random: maximum error: 
OK\? : YES
After matrix solving: maximum error: 
OK\? : YES
EOT


#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("circ/check_circ.sal","Constraint bug, cant have null rows in idmap class",<<'EOT');
error in u symmetry = .* OK . 1
error in v symmetry = .* OK . 1
error in h symmetry = .* OK . 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("bug100/save.state.tmp",
       "Advdif check bug100 is fixed",<<'EOT');
1.230.*e-01
1.230.*e-01
__NO_SKIP__
1.249.*e-01
1.249.*e-01
1.230.*e-01
1.230.*e-01
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("newff/newff.check_fields.out",
       "Check fields<ndof in constraints",<<'EOT');
Assertion failed: 
read_mesh: Read field
newff.depl:
PETSC-FEM error 
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("newff/newff.check_nodes.out",
       "Check node<nnod in constraints",<<'EOT');
Assertion failed: "node<=nnod."
read_mesh: Read node= 35 greater than nnod= 34
newff.depl:.*: "1. 35 1 \-1. 33 1"
---------------
PETSC-FEM error at file "readmesh.cpp", line 
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("lupart/check_part.verif.tmp",
       "Partitioning test",<<'EOT');
Random partitioning OK \? > 1
EOT


#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/distmaps.sal.tmp",
       "Distributed map<int,double>, sched alg NOT GROUPED",<<'EOT');
Args: .* sched 0
error < tol OK \? > 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/distcont.sal.tmp",
       "Distributed container class.",read_file("distcont.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/distmapg.sal.tmp",
       "Distributed map<int,double>, sched alg GROUPED",<<'EOT');
Args: .* sched 1
error < tol OK \? > 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/distmat.sal.tmp",
       "Distributed matrix<int,int,double>",<<'EOT');
error < tol OK \? yes
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("lupart/check_direct_superlu.verif.tmp",
       "SuperLU direct solver (SparseDirect class)",<<'EOT');
Direct/SuperLU  OK \? > 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("lupart/check_direct_petsc.verif.tmp",
       "PETSc direct solver (SparseDirect class)",<<'EOT');
Direct/PETSc  OK \? > 1, 
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sqcav/output.CASE_sqcav.np_1.case_lu.out.tmp",
       "Measure performance test.",read_file("test_meas_perf.ans.txt"));


#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("nsgenload/verif.constraint.tmp",
       "Cyclic reference constraints ",<<'EOT');
Error < tol OK \? 1
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Laplace tests. ');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sector/sector_chunk1.sal",
       "Periodic case with chunk_size = 1 ",$sector_test);
expect("sector/sector_chunk10.sal",
       "Periodic case with chunk_size = 10 ",$sector_test);
expect("sector/sector_chunk100.sal",
       "Periodic case with chunk_size = 100 ",$sector_test);

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sector/sector_triangle_npg3.sal",
       "Periodic case with triangles npg=3 ",$sector_test);
expect("sector/sector_triangle_npg4.sal",
       "Periodic case with triangles npg=4 ",$sector_test);
expect("sector/sector_triangle_npg7.sal",
       "Periodic case with triangles npg=7 ",$sector_test);

end_section();

begin_section('Misc tests.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fstack/tfstack.out","FileStack class",<<'EOT');
line 1 on file_stack 1: <1 line continued ...    continued ....>
replaced line 2 on file_stack 1: <replaced line>
line 3 on file_stack 2: <1 line continued ...    continued .... .copy 2.>
line 4 on file_stack 1: <2 several blank lines ...>
replaced line 5 on file_stack 1: <replaced line>
line 6 on file_stack 2: <2 several blank lines ... .copy 2.>
line 7 on file_stack 1: <1b line continued ...    continued ....>
replaced line 8 on file_stack 1: <replaced line>
line 9 on file_stack 2: <1b line continued ...    continued .... .copy 2.>
line 10 on file_stack 1: <2b several blank lines ...>
replaced line 11 on file_stack 1: <replaced line>
line 12 on file_stack 2: <2b several blank lines ... .copy 2.>
line 13 on file_stack 1: <3b third line>
replaced line 14 on file_stack 1: <replaced line>
line 15 on file_stack 2: <3b third line .copy 2.>
line 16 on file_stack 1: <1a line continued ...    continued ....>
replaced line 17 on file_stack 1: <replaced line>
line 18 on file_stack 2: <1a line continued ...    continued .... .copy 2.>
line 19 on file_stack 1: <2a several blank lines ...>
replaced line 20 on file_stack 1: <replaced line>
line 21 on file_stack 2: <2a several blank lines ... .copy 2.>
line 22 on file_stack 1: <3a third line>
replaced line 23 on file_stack 1: <replaced line>
line 24 on file_stack 2: <3a third line .copy 2.>
line 25 on file_stack 1: <3 third line>
replaced line 26 on file_stack 1: <replaced line>
line 27 on file_stack 2: <3 third line .copy 2.>
cat file1 reversed at the start of file2:
line 1 : <3 third line .copy 2.>
line 2 : <3a third line .copy 2.>
line 3 : <2a several blank lines ... .copy 2.>
line 4 : <1a line continued ...    continued .... .copy 2.>
line 5 : <3b third line .copy 2.>
line 6 : <2b several blank lines ... .copy 2.>
line 7 : <1b line continued ...    continued .... .copy 2.>
line 8 : <2 several blank lines ... .copy 2.>
line 9 : <1 line continued ...    continued .... .copy 2.>
line 10 : <1 line continued ...    continued ....>
line 11 : <2 several blank lines ...>
line 12 : <1b line continued ...    continued ....>
line 13 : <2b several blank lines ...>
line 14 : <3b third line>
line 15 : <1a line continued ...    continued ....>
line 16 : <2a several blank lines ...>
line 17 : <3a third line>
line 18 : <3 third line>
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fstack/fstack2.verif.tmp","FileStack class, file name and line pos.",<<'EOT');
file ok. 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("texthash/tthash.sal","TextHashTable class",<<'EOT');
Global counting OK
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("texthash/thash2.output.tmp","TextHashTable::find",
                     read_file('texthash/thash2.ans'));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("genload/output.case_fstack.tmp","Report error for bad number of nodes in conn. line",<<'EOT');
Assertion failed:
Error reading .* at
genload.depl:.*"35 36 37 38"
PETSC-FEM error at file "readmesh.cpp"
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sparse/output.sparse.superlu.tmp","Sparse Mat/Vec classes, SuperLU version",
                        read_file("sparse.test"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("sparse/output.sparse.petsc.tmp","Sparse Mat/Vec classes, PETSc version",
                        read_file("sparse.test"));

end_section();

begin_section('Advdif tests');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advdif/sine_fine_mesh.out","Advdif // conv. to analytic in fine mesh",<<'EOT');
Dt=.* error=.* \< tol_error=.* OK\?.*1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advdif/sine_crank_nic.out","Advdif // quad. conv. for Crank Nic.",<<'EOT');
||u_16-u_128|| = 
||u_32-u_128|| = 
||u_32-u_128|| / ||u_16-u_128|| = .*, < 0.25 OK\? 1 
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Turbchan tests.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("turbchan/turbchan.ver.tmp","Advdif/Shallw Water Turb. // flow in a channel",<<'EOT');
Residual at iteration 10 =.*OK : 1
symmetry OK : 1
transversal velocity small OK : 1
k at output .* Correct value OK : 1
e at output .* Correct value OK : 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("turbchanw/heat.out.tmp","Advdif/Shallw Water Turb. // log vars // heat eq.",<<'EOT');
2.397
2.312
2.303
2.302
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("turbchanw/check.swturb.verif.tmp",
        "Advdif/Shallw Water Turb. // log vars // shallow w. eq.",<<'EOT');
Max rel error OK
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('New flux function tests.');

sub check_newff {
    my $case=shift();
    expect("newff/$case.vrf.tmp","newff/test case $case",
	   'rel. max error  <tol\(.*\) OK\? 1');
}

check_newff('adv_dif_stdy_bcconv');
check_newff('adv_dif_temp');
check_newff('adv_temp_cp');
check_newff('dif_temp');
check_newff('dif_temp_cp');
check_newff('full_full_jacs');
check_newff('full_full_jacs_cp');
check_newff('full_full_jacs_adv_cp');
check_newff('full_full_jacs_t');
check_newff('full_full_jacs_tr');
check_newff('full_full_jacs_wf');
check_newff('full_jacs');
check_newff('full_jacs_cp');
check_newff('pure_adv');
check_newff('reac_adv_dif_temp_y');
check_newff('reac_dif_temp');
check_newff('reac_steady');
check_newff('std_ard_x_y');
check_newff('stdy_dif');

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Burgers/advdif tests.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("burgers/save.state._d01_wf.tmp",
       "Advdif/Burgers with weak form+bcconv, d=0.1",<<'EOT');
9.64.*e\-01  
7.59.*e\-01  
e\-15
-7.59.*e\-01  
-9.64.*e\-01  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("burgers/save.state._d01_nwf.tmp",
       "Advdif/Burgers with no weak form, d=0.1",<<'EOT');
9.64.*e\-01  
7.59.*e\-01  
e\-15
-7.59.*e\-01  
-9.64.*e\-01  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("burgers/test__d01.out.tmp",
       "Advdif/Burgers coincidence weak/no weak form, d=0.1",<<'EOT');
max. error .*, <tol\(.*\) OK\? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("burgers/save.state._d001_wf.tmp",
       "Advdif/Burgers with weak form+bcconv, d=0.01",<<'EOT');
1.000.*e\+00  
1.000.*e\+00  
1.094.*e\+00  
__NO_SKIP__
1.094.*e\+00  
e\-1(4|5)  
e\-1(4|5)  
-1.094.*e\+00  
-1.094.*e\+00  
-1.000.*e\+00  
-1.000.*e\+00  
EOT

=cut
1.00.*e\+00  
1.00.*e\+00  
1.09.*e\+00  
__NO_SKIP__
1.09.*e\+00  
e\-1(4|5|6)
e\-1(4|5|6)
-1.09.*e\+00  
-1.09.*e\+00  
-1.00.*e\+00  
-1.00.*e\+00  
=cut

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("burgers/save.state._d001_nwf.tmp",
       "Advdif/Burgers with no weak form, d=0.01",<<'EOT');
1.000.*e\+00  
1.000.*e\+00  
1.09(3|4).*e\+00  
__NO_SKIP__
1.09(3|4).*e\+00  
e\-1(4|5|6|7)
e\-1(4|5|6|7)
-1.09(3|4).*e\+00  
-1.09(3|4).*e\+00  
-1.000.*e\+00  
-1.000.*e\+00  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("burgers/test__d001.out.tmp",
       "Advdif/Burgers coincidence weak/no weak form, d=0.01",<<'EOT');
max. error .*, <tol\(.*\) OK\? 1
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Linear advection (advdif) tests.');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/skew.verif.tmp",
       "Steady advection skew to the mesh",<<'EOT');
Over shoot < tol OK \? 1
Under shoot < tol OK \? 1
Number of nodes with error > tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/parallel_cone.verif.tmp",
       "Unsteady advection of cone parallel to mesh",<<'EOT');
Error in x position < tol OK \? 1
Error in y position < tol OK \? 1
Max grouth in x second moment OK \? 1
Max grouth in y second moment OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/skew_cone.verif.tmp",
       "Unsteady advection of cone skew to mesh",<<'EOT');
Error in x position < tol OK \? 1
Error in y position < tol OK \? 1
Max grouth in x second moment OK \? 1
Max grouth in y second moment OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/rotating_cone.verif.tmp",
       "Unsteady advection of rotating cone",<<'EOT');
Error in x position < tol OK \? 1
Error in y position < tol OK \? 1
Max grouth in x second moment OK \? 1
Max grouth in y second moment OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/skew_noise.verif.tmp",
       "Steady advection skew to the mesh",<<'EOT');
Over shoot < tol OK \? 1
Under shoot < tol OK \? 1
Number of nodes with error > tol OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/parallel_cone_noise.verif.tmp",
       "Unsteady advection of cone parallel to mesh",<<'EOT');
Error in x position < tol OK \? 1
Error in y position < tol OK \? 1
Max grouth in x second moment OK \? 1
Max grouth in y second moment OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/skew_cone_noise.verif.tmp",
       "Unsteady advection of cone skew to mesh",<<'EOT');
Error in x position < tol OK \? 1
Error in y position < tol OK \? 1
Max grouth in x second moment OK \? 1
Max grouth in y second moment OK \? 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("advec/rotating_cone_noise.verif.tmp",
       "Unsteady advection of rotating cone",<<'EOT');
Error in x position < tol OK \? 1
Error in y position < tol OK \? 1
Max grouth in x second moment OK \? 1
Max grouth in y second moment OK \? 1
EOT

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Generic Load element');

$genl_check = <<'EOT';
5.62.*e\-01  5.562.*e\+00  
__NO_SKIP__
5.62.*e-01  5.562.*e\+00  
.
.
.
.
-5.62.*e-01  4.43.*e\+00  
-5.62.*e-01  4.43.*e\+00  
__SKIP__
-8.7.*e-01  4.12.*e\+00  
__NO_SKIP__
-8.7.*e-01  4.12.*e\+00  
-9.37.*e-01  4.062.*e\+00  
-9.37.*e-01  4.062.*e\+00  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("genload/save.case_1.tmp",
       "Generic Load element one gen.load.- double".
       " layer hfilm elem. at the end",$genl_check);
expect("genload/save.case_2a.tmp",
       "Generic Load element one gen.load.- double".
       " layer with source elem. at the end",$genl_check);
expect("genload/save.case_2b.tmp",
       "Generic Load element one gen.load.- single".
       " layer with source elem. at the end",$genl_check);
expect("genload/save.case_genl1d.tmp",
       "Generic Load, 0d element with source term",$genl_check);
expect("genload/save.case_genl1dh.tmp",
       "Generic Load, 0d element with hfilm coeff",$genl_check);

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('IISD solver');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("lupart/check_iisd.verif.tmp",
       "IISD solver",<<'EOT');
IISD on 1 processors OK \? > 1, 
IISD on 2 processors OK \? > 1, 
IISD on 2 processors with rand part. OK \? > 1
IISD on 2 processors with CGS OK \? > 1, 
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#  expect("lupart/check_iisd_superlu.verif.tmp",
#         "IISD/SuperLU solver",<<'EOT');
#  __EXACT_MATCH__
#  IISD/SuperLU on 1 processors OK ? > 1, 
#  IISD/SuperLU on 2 processors OK ? > 1, 
#  IISD/SuperLU on 2 processors with rand part. OK ? > 1, 
#  IISD/SuperLU on 2 processors with CGS OK ? > 1
#  EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("lupart/check_peri.verif.tmp",
       "IISD solver",<<'EOT');
IISD on 2 processors with periodic b.c.'s OK \? > 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("lupart/check_subpart.verif.tmp",
       "IISD/SubPartitioning solver",<<'EOT');
IISD/Subpartitioning  OK \? > 1
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
for ($case=1; $case<=4; $case++) {
#    print "checking pfmat/output.case${case}iip.tmp\n";
    expect("pfmat/output.case${case}iip.tmp",
          "PFMat/case$case/(local_solver=petsc)",
          "All tests OK.*1");
}

=cut
# Unfortunately, this doesn't work.
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
for ($case=1; $case<=4; $case++) {
#    print "checking pfmat/output.case${case}iip.tmp\n";
    expect("pfmat/output.case${case}iip.tmp",
          "PFMat/case$case/(local_solver=petsc)",
          "All tests OK.*1");
}
=cut

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
for ($case=1; $case<=2; $case++) {
#    print "checking pfmat/output.case${case}p.tmp\n";
    expect("pfmat/output.case${case}iip.tmp",
          "PETScMat/case$case",
          "All tests OK.*1");
}

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("pfmat/output.case1sdp.tmp",
       "SparseDirect(solver PETSc)",
       "All tests OK.*1");

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("pfmat/output.case1sds.tmp",
       "SparseDirect(solver SuperLU)",
       "All tests OK.*1");

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('IISD solver. Checking with advdif');

check_newff('adv_dif_temp_iisd');
check_newff('adv_temp_cp_iisd');
check_newff('dif_temp_iisd');
check_newff('full_full_jacs_iisd');
check_newff('full_jacs_iisd');
check_newff('full_jacs_cp_iisd');
check_newff('pure_adv_iisd');
check_newff('reac_adv_dif_temp_y_iisd');
check_newff('std_ard_x_y_iisd');

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Graph Partitioning');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/dist_graph.np1.out.tmp","Distributed graph 1proc.",
               read_file("distmap/dist_graph.np1.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/dist_graph.np2.out.tmp","Distributed graph 2proc.",
               read_file("distmap/dist_graph.np2.ans.txt"));

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("distmap/dist_graph.np3.out.tmp","Distributed graph 3proc.",
               read_file("distmap/dist_graph.np3.ans.txt"));

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Advective/Hydrology module');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/plain.verif.tmp","Flat aquifer bottom",
       'test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/var_eta0.verif.tmp","Flat aquifer bottom as varying elem prop.",
       'test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/const.verif.tmp","Flat aquifer bottom coincide.",
       'OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/var_eta1.verif.tmp","Variable aquifer bottom",
       'test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/transient.verif.tmp","Transient - Constant aquifer bottom",
       'test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/stream.verif.tmp","Stream elemset",'Test OK. > 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/aquist.verif.tmp","aquifer+stream coupling",'Test OK. > 1');

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

begin_section('Dynamically loadable amplitude functions');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/output.dl_fun1.verif.tmp","Linear ramp Dirichlet b.c.",
       'Test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/output.dl_fun2.verif.tmp","Smooth ramp Dirichlet b.c.",
       'Test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/output.dl_fun3.verif.tmp","Tanh ramp Dirichlet b.c.",
       'Test OK \? 1');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("aquifer/var_rain.verif.tmp","Variable test rain.",
       'Variable rain test OK \? 1');

end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 

print "\n",'-' x 50,"\n\n";
end_section();

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
/'/; # to disable auto-filling an indenting in Emacs

final_check();
