global N x3
source("data.m.tmp");
u = aload("cubcav.iisd_sbp1_uj.tmp");

x3=aload("cubcav.nod.tmp");
rows(x)==(N+1)^3 || error("Don't match state length");

x=x3(1:N+1,1);
y=x;
z=x;

rem(N,2)==0 || error("N should be pair");
N2 = N/2;

