// This file was automatically created by "mkcode.pl"
// It is supposed not to be modified by hand. 
// It is not a standard header file. You should
// not include it in a program, this is only
// included once in file "fm2prod2.cpp"
//--------------------------------------------------------------------------------
DEFFUN2(p_1_1_1) {
c[0] = a[0]*b[0];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_1_2) {
c[0] = a[0]*b[0];
c[1] = a[0]*b[1];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_1_3) {
c[0] = a[0]*b[0];
c[1] = a[0]*b[1];
c[2] = a[0]*b[2];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_2_1) {
c[0] = a[0]*b[0]+a[1]*b[1];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_2_2) {
c[0] = a[0]*b[0]+a[1]*b[2];
c[1] = a[0]*b[1]+a[1]*b[3];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_2_3) {
c[0] = a[0]*b[0]+a[1]*b[3];
c[1] = a[0]*b[1]+a[1]*b[4];
c[2] = a[0]*b[2]+a[1]*b[5];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_3_1) {
c[0] = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_3_2) {
c[0] = a[0]*b[0]+a[1]*b[2]+a[2]*b[4];
c[1] = a[0]*b[1]+a[1]*b[3]+a[2]*b[5];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_1_3_3) {
c[0] = a[0]*b[0]+a[1]*b[3]+a[2]*b[6];
c[1] = a[0]*b[1]+a[1]*b[4]+a[2]*b[7];
c[2] = a[0]*b[2]+a[1]*b[5]+a[2]*b[8];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_1_1) {
c[0] = a[0]*b[0];
c[1] = a[1]*b[0];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_1_2) {
c[0] = a[0]*b[0];
c[1] = a[0]*b[1];
c[2] = a[1]*b[0];
c[3] = a[1]*b[1];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_1_3) {
c[0] = a[0]*b[0];
c[1] = a[0]*b[1];
c[2] = a[0]*b[2];
c[3] = a[1]*b[0];
c[4] = a[1]*b[1];
c[5] = a[1]*b[2];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_2_1) {
c[0] = a[0]*b[0]+a[1]*b[1];
c[1] = a[2]*b[0]+a[3]*b[1];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_2_2) {
c[0] = a[0]*b[0]+a[1]*b[2];
c[1] = a[0]*b[1]+a[1]*b[3];
c[2] = a[2]*b[0]+a[3]*b[2];
c[3] = a[2]*b[1]+a[3]*b[3];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_2_3) {
c[0] = a[0]*b[0]+a[1]*b[3];
c[1] = a[0]*b[1]+a[1]*b[4];
c[2] = a[0]*b[2]+a[1]*b[5];
c[3] = a[2]*b[0]+a[3]*b[3];
c[4] = a[2]*b[1]+a[3]*b[4];
c[5] = a[2]*b[2]+a[3]*b[5];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_3_1) {
c[0] = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
c[1] = a[3]*b[0]+a[4]*b[1]+a[5]*b[2];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_3_2) {
c[0] = a[0]*b[0]+a[1]*b[2]+a[2]*b[4];
c[1] = a[0]*b[1]+a[1]*b[3]+a[2]*b[5];
c[2] = a[3]*b[0]+a[4]*b[2]+a[5]*b[4];
c[3] = a[3]*b[1]+a[4]*b[3]+a[5]*b[5];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_2_3_3) {
c[0] = a[0]*b[0]+a[1]*b[3]+a[2]*b[6];
c[1] = a[0]*b[1]+a[1]*b[4]+a[2]*b[7];
c[2] = a[0]*b[2]+a[1]*b[5]+a[2]*b[8];
c[3] = a[3]*b[0]+a[4]*b[3]+a[5]*b[6];
c[4] = a[3]*b[1]+a[4]*b[4]+a[5]*b[7];
c[5] = a[3]*b[2]+a[4]*b[5]+a[5]*b[8];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_1_1) {
c[0] = a[0]*b[0];
c[1] = a[1]*b[0];
c[2] = a[2]*b[0];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_1_2) {
c[0] = a[0]*b[0];
c[1] = a[0]*b[1];
c[2] = a[1]*b[0];
c[3] = a[1]*b[1];
c[4] = a[2]*b[0];
c[5] = a[2]*b[1];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_1_3) {
c[0] = a[0]*b[0];
c[1] = a[0]*b[1];
c[2] = a[0]*b[2];
c[3] = a[1]*b[0];
c[4] = a[1]*b[1];
c[5] = a[1]*b[2];
c[6] = a[2]*b[0];
c[7] = a[2]*b[1];
c[8] = a[2]*b[2];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_2_1) {
c[0] = a[0]*b[0]+a[1]*b[1];
c[1] = a[2]*b[0]+a[3]*b[1];
c[2] = a[4]*b[0]+a[5]*b[1];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_2_2) {
c[0] = a[0]*b[0]+a[1]*b[2];
c[1] = a[0]*b[1]+a[1]*b[3];
c[2] = a[2]*b[0]+a[3]*b[2];
c[3] = a[2]*b[1]+a[3]*b[3];
c[4] = a[4]*b[0]+a[5]*b[2];
c[5] = a[4]*b[1]+a[5]*b[3];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_2_3) {
c[0] = a[0]*b[0]+a[1]*b[3];
c[1] = a[0]*b[1]+a[1]*b[4];
c[2] = a[0]*b[2]+a[1]*b[5];
c[3] = a[2]*b[0]+a[3]*b[3];
c[4] = a[2]*b[1]+a[3]*b[4];
c[5] = a[2]*b[2]+a[3]*b[5];
c[6] = a[4]*b[0]+a[5]*b[3];
c[7] = a[4]*b[1]+a[5]*b[4];
c[8] = a[4]*b[2]+a[5]*b[5];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_3_1) {
c[0] = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
c[1] = a[3]*b[0]+a[4]*b[1]+a[5]*b[2];
c[2] = a[6]*b[0]+a[7]*b[1]+a[8]*b[2];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_3_2) {
c[0] = a[0]*b[0]+a[1]*b[2]+a[2]*b[4];
c[1] = a[0]*b[1]+a[1]*b[3]+a[2]*b[5];
c[2] = a[3]*b[0]+a[4]*b[2]+a[5]*b[4];
c[3] = a[3]*b[1]+a[4]*b[3]+a[5]*b[5];
c[4] = a[6]*b[0]+a[7]*b[2]+a[8]*b[4];
c[5] = a[6]*b[1]+a[7]*b[3]+a[8]*b[5];
}
//--------------------------------------------------------------------------------
DEFFUN2(p_3_3_3) {
c[0] = a[0]*b[0]+a[1]*b[3]+a[2]*b[6];
c[1] = a[0]*b[1]+a[1]*b[4]+a[2]*b[7];
c[2] = a[0]*b[2]+a[1]*b[5]+a[2]*b[8];
c[3] = a[3]*b[0]+a[4]*b[3]+a[5]*b[6];
c[4] = a[3]*b[1]+a[4]*b[4]+a[5]*b[7];
c[5] = a[3]*b[2]+a[4]*b[5]+a[5]*b[8];
c[6] = a[6]*b[0]+a[7]*b[3]+a[8]*b[6];
c[7] = a[6]*b[1]+a[7]*b[4]+a[8]*b[7];
c[8] = a[6]*b[2]+a[7]*b[5]+a[8]*b[8];
}
