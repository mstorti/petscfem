##__INSERT_LICENSE__
## $Id: loadiisd.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
petsc_data_name="system.dat";
petscload
resiisd=res;
dxiisd=dx;

#petsc_data_name="debug_iisd.dat";
#petscload

load -force map.dat

s = glob("a_ll_*");

for ss = s'
  ss=ss';
  petsc_data_name=ss;
  petscload
  eval([ss "=getblock(" ss ");"]);
endfor
nproc = rows(s);

A_LI=getblock(atet);
A_IL=getblock(atet_0);
A_II=getblock(atet_1);

clear atet atet_0 atet_1

map=map(:,2);
