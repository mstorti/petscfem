## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");

w=zhomo([0 h 0 L],2,N+1,[1 0 1 1 xratio 1]);
[x2,i2]=pfcm2fem(w);
i2 = i2(:,[1 4 3 2]);

nnod = rows(x2);		# number of real nodes
nelem = rows(i2);

fid = fopen("strip2d.bc.tmp","w");

NN = N+1;			# number of nodes in one dimension

## impose pressure in some real node
fprintf(fid,"%d %d 0.\n",1,3);
## No-slip condition at both ends
for node=[1 NN]
  for dof=1:3
    fprintf(fid,"%d %d 0.\n",node,dof);
  endfor
endfor

fid = fopen("strip2d.cnstr.tmp","w");
## Periodic b.c.'s 
for k=1:NN
  for dof=1:3
    fprintf(fid,"-1. %d %d 1. %d %d\n",k+NN,dof,k,dof);
  endfor
endfor

fclose(fid);

asave("strip2d.nod.tmp",x2);
asave("strip2d.con.tmp",i2);
