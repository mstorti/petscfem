## __INSERT_LICENSE__
## $Id: mkstrip.m,v 1.2 2003/01/10 19:17:09 mstorti Exp $

source("data.m.tmp");

w=zhomo([0 h 0 L],2,N+1,[1 0 1 1 xratio 1]);
[x2,i2]=pfcm2fem(w);
i2 = i2(:,[1 4 3 2]);
[x3,i3]=extrude(x2,i2,1,h);

nnod = rows(x3);		# number of real nodes
nodf_u = nnod+1;                # fictitious node for ux,uy
nodf_r = nnod+2;                # fictitious node for rotations
nfic = 2;			# number of fictitious nodes

asave("strip.fic.tmp",[nodf_u; nodf_r]);

x3=[x3;
    0 0 0;
    0 0 0];

nelem = rows(i3);
i3 = [i3, nodf_u*ones(nelem,1), nodf_r*ones(nelem,1)];

fid = fopen("strip.bc.tmp","w");

NN = N+1;			# number of nodes in one dimension

## impose pressure in some real node
fprintf(fid,"%d %d 0.\n",1,4);
## No-slip condition at both ends
for node=[1 NN]
  for dof=1:3
    fprintf(fid,"%d %d 0.\n",node,dof);
  endfor
endfor

## Impose k-epsilon to some value
tke=0.1;
for k=1:NN
  for dof=5:6
    fprintf(fid,"%d %d %f\n",k,dof,tke);
  endfor
endfor
fclose(fid);

fid = fopen("strip.cnstr.tmp","w");
## Periodic b.c.'s 
for k=1:NN
  for dof=1:6
    fprintf(fid,"-1. %d %d 1. %d %d\n",k+NN,dof,k,dof);
    fprintf(fid,"-1. %d %d 1. %d %d\n",k+2*NN,dof,k,dof);
    fprintf(fid,"-1. %d %d 1. %d %d\n",k+3*NN,dof,k,dof);
  endfor
endfor

fclose(fid);

asave("strip.nod.tmp",x3);
asave("strip.con.tmp",i3);
i3p = [i3(:,[1 5 6 4 8 7]);
       i3(:,[6 2 1 7 3 4])];
i3p = [i3p ones(rows(i3p),2)*diag([nodf_u nodf_r])];
asave("strip.con-prism.tmp",i3p);
