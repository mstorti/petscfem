#!/usr/bin/perl

require 'myexpect.pl';

##------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
/'/; # to disable auto-filling an indenting in Emacs
expect("fastmat2a.sal","FastMat2 library",<<'EOT');
--- A:  ----
 11 12 13
 21 22 23
 31 32 33
free indices  3
--- A first row
 11 12 13
--- A:  ----
 11 12 13
 21 22 23
 31 32 33
--- A column 2:  ----
 12 22 32
--- Column 3 ----
 13
 23
 33
--- all the matrix:  ----
 11 12 13
 21 22 23
 31 32 33
--- Column 3 \(but as row now\) ----
 13 23 33
--- all the matrix:  ----
 11 12 13
 21 22 23
 31 32 33
--- A:  ----
 1 2 0 0 0
 3 4 0 0 0
 0 0 0 0 0
 0 0 0 0 0
 0 0 0 0 0
--- Adds to elements
 2 4 0 0 0
 6 8 0 0 0
 0 0 0 0 0
 0 0 0 0 0
 0 0 0 0 0
--- Scales
 1 2 0 0 0
 3 4 0 0 0
 0 0 0 0 0
 0 0 0 0 0
 0 0 0 0 0
--- Sub-block
 4 3
 2 1
--- A:  ----
 4 2
 3 1
--- B:  ----
 4 2
 3 1
--- B \(reshaped\):  ----
 4 2 3 1
--- A':  ----
 4 3
 2 1
--- B':  ----
 4 3
 2 1
B
 8 6
 4 2
B
 4 3
 2 1
B
 44 33
 22 11
--- 4
 12 10
 8 6
A set to 2
 2
 2
 2
 2
 2
--- A:  ----
 1 1 1
 1 1 1
--- B:  ----
 2 2
 2 2
 2 2
--- C = A\*B 
 6 6
 6 6
--- Kron prod 
 2 2
 2 2
 2 2
--- D:  ----
for first indices ->   1  1  
 4 4 4
 4 4 4
for first indices ->   1  2  
 4 4 4
 4 4 4
for first indices ->   2  1  
 4 4 4
 4 4 4
for first indices ->   2  2  
 4 4 4
 4 4 4
for first indices ->   3  1  
 4 4 4
 4 4 4
for first indices ->   3  2  
 4 4 4
 4 4 4
--- D_{ijki}:  ----
 12 12
 12 12
--- F_{lm} = D_{ijkl} D_{ijkm} ----
 192 192 192
 192 192 192
 192 192 192
--- E = D_{ijjk} =  ----
 12 12
 12 12
--- H = Sum over first index of G:  ----
 -6 -6
 -6 -6
--- H = Sum_square over first index of G:  ----
 18 18
 18 18
--- H = Sum_abs over columns of G:  ----
 6 6
 6 6
--- H = Min over first index of G:  ----
 -5 -3
 -3 -3
--- H = Max over first index of G:  ----
 -3 -3
 -3 5
--- H = Min abs over first index of G:  ----
 3 3
 3 3
--- H = Max abs over first index of G:  ----
 5 3
 3 5
diag\(K\)
 111 122 133
 211 222 233
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fastmat2.sal","FastMat2 library with caches",<<'EOT');
--- A:  ----
 1 2 13
 3 4 23
 31 32 33
--- B:  ----
 11 24 169
 63 88 529
 961 1024 1089
--- C:  ----
 1.5 2.5 13
 3.5 4.5 23
 31 32 33
--- D:  ----
 0.5 1 13
 1.5 2 23
 31 32 33
--- E:  ----
 1 2 4
 5 3 2
 3 2 1
--- G:  ----
 16 30 96
--- K:  ----
 1 2 13
 3 -1000 23
 31 32 33
--- L:  ----
 16 1026 96
--- M:  ----
 31 32 33
--- P:  ----
 1 -1000 13
--- Q:  ----
 31 1000 33
--- R:  ----
 1 2 13
--- S:  ----
 3.14 3.14 3.14
 3.14 3.14 3.14
 3.14 3.14 3.14
--- U:  ----
 2 4 26
 6 8 46
 62 64 66
--- V:  ----
 6 7 18
 8 9 28
 36 37 38
--- W:  ----
 410 426 488
 728 758 890
 1150 1246 2228
--- X:  ----
 11 12 13 21 22 23 31 32 33
--- Y:  ----
 1 3 5
 7 9 11
 13 15 17
NA: 1.000000 3.000000 5.000000 
7.000000 9.000000 11.000000 
13.000000 15.000000 17.000000 

NB: 1.000000 3.000000 5.000000 
7.000000 9.000000 11.000000 
13.000000 15.000000 17.000000 

--- Z:  ----
not defined!
--- Z1:  ----
 1 4 33
--- Z2:  ----
 1 3 2
--- Z3:  ----
38.000000
double(Z3): 38
--- Z5:  ----
 13 14
 23 24
 33 34
 43 44
 53 54
 23 24
 33 34
 43 44
 53 54
 63 64
 73 74
 83 84
 93 94
 103 104
 113 114
 123 124
 133 134
 143 144
 153 154
 163 164
 173 174
 183 184
 193 194
 203 204
--- Z6:  ----
not defined!
--- Z8:  ----
 1 2 4
 5 3 2
 3 2 1
dz8: 5
--- Z9:  ----
 -0.2 1.2 -1.6
 0.2 -2.2 3.6
 0.2 0.8 -1.4
--- Z11:  ----
 1 2 2 4 13 26
 3 4 6 8 39 52
 3 6 4 8 23 46
 9 12 12 16 69 92
 31 62 32 64 33 66
 93 124 96 128 99 132
--- Z12:  ----
not defined!
z12: 3802
z121: 3802
z13: 33
z14: 1
z15: 180
--- Z16:  ----
 0.111111 -5.68651e-17 -0.333333 0.222222
 -1.81641e-17 -0.5 1 -0.5
 -0.333333 1 -5.4 3.73333
 0.222222 -0.5 4.63333 -3.35556
--- Z17:  ----
 1.11752 1 0.716531 1.24885
 1 0.606531 2.71828 0.606531
 0.716531 2.71828 0.00451658 41.8183
 1.24885 0.606531 102.856 0.03489
--- Z18:  ----
 0.000123457 1.98271e-34 0.00111111 0.000493827
 4.53159e-37 0.0025 0.01 0.0025
 0.00111111 0.01 0.2916 0.139378
 0.000493827 0.0025 0.214678 0.112598
z115: 180
--- Z116:  ----
 0.111111 -5.68651e-17 -0.333333 0.222222
 -1.81641e-17 -0.5 1 -0.5
 -0.333333 1 -5.4 3.73333
 0.222222 -0.5 4.63333 -3.35556
--- Z117:  ----
 1.11752 1 0.716531 1.24885
 1 0.606531 2.71828 0.606531
 0.716531 2.71828 0.00451658 41.8183
 1.24885 0.606531 102.856 0.03489
--- Z118:  ----
 0.000123457 3.23364e-35 0.00111111 0.000493827
 3.29935e-36 0.0025 0.01 0.0025
 0.00111111 0.01 0.2916 0.139378
 0.000493827 0.0025 0.214678 0.112598
--- Z19:  ----
 3.342 0 0
 0 3.342 0
 0 0 3.342
--- Z21:  ----
 11 12 13 14 15
 21 22 23 24 25
 31 32 33 34 35
--- Z21:  ----
 11 12 13 14 15
 21 22 23 24 25
 31 32 33 34 35
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fastmat2c.sal","FastMat2 library. FEM problem.",<<'EOT');
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

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
# Here we check only the $y$ components of velocity, since the others are
# rather spurious. 
expect("oscplate1.sal","Time dependent boundary conditions",<<'EOT');
5.100.*e-01
1.819.*e-01
1.980.*e-02
-3.595.*e-02
-4.028.*e-02
-2.769.*e-02
-1.451.*e-02
-5.920.*e-03
-1.748.*e-03
EOT

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
norm\(u16-u128\): 0.011
norm\(u32-u128\): 0.002
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

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("tfastvec.sal",
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


#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
expect("fstack/tfstack.sal","FileStack class",<<'EOT');
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
expect("texthash/tthash.sal","TextHashTable class",<<'EOT');
Global counting OK
EOT

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
expect("turbchanw/swturb.out.tmp","Advdif/Shallw Water Turb. // log vars // shallow w. eq.",<<'EOT');
3\.13.*e-02 .*  1\.00.*e-01  -7\.95.*e\+00  -9\.50.*e\+00  
3\.37.*e-02 .*  1\.00.*e-01  -8\.02.*e\+00  -9\.62.*e\+00  
3\.40.*e-02 .*  1\.00.*e-01  -8\.17.*e\+00  -9\.86.*e\+00  
3\.42.*e-02 .*  1\.00.*e-01  -8\.29.*e\+00  -1\.00.*e\+01  
3\.44.*e-02 .*  1\.00.*e-01  -8\.39.*e\+00  -1\.01.*e\+01  
3\.46.*e-02 .*  1\.00.*e-01  -8\.46.*e\+00  -1\.02.*e\+01  
3\.46.*e-02 .*  1\.00.*e-01  -8\.50.*e\+00  -1\.02.*e\+01  
3\.47.*e-02 .*  1\.00.*e-01  -8\.53.*e\+00  -1\.02.*e\+01  
EOT

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
sub check_newff {
    my $case=shift();
    expect("newff/$case.vrf.tmp","newff/test case $case",
	   'rel. max error  <tol\(.*\) OK\? 1');
}

check_newff('dif_temp');
check_newff('adv_dif_temp');
check_newff('reac_adv_dif_temp_y');
check_newff('reac_steady');
check_newff('reac_dif_temp');
check_newff('std_ard_x_y');
check_newff('full_jacs');
check_newff('full_full_jacs');
check_newff('full_full_jacs_wf');

#------/*/------/*/------/*/------/*/------/*/------/*/------/*/ 
/'/; # to disable auto-filling an indenting in Emacs

print "\n",'-' x 80,"\n\n";

final_check();
