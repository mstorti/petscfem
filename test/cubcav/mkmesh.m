##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.9 2004/01/30 02:24:02 mstorti Exp $
source("data.m.tmp");
#  N=5;
#  hratio = 0;
#  tol=1e-3/N;
#  leaky = 0;

w=zhomo([0 1 0 1],N+1,N+1,[1 hratio 1 1 hratio 1]);
[x2,i2] = pfcm2fem(w);
x2 = x2(:,[2 1]);
xx=x2(1:N+1,1);

[x3,i3] = extrude(x2,i2,N,xx);

x = x3(:,1);
y = x3(:,2);
z = x3(:,3);

wall_f = abs(x)<tol  | abs(1-x)<tol | abs(z)<tol | abs(1-z)<tol | \
     abs(y)<tol;

top_f = abs(y-1)<tol;

if leaky 
  wall = find(wall_f & ! top_f);
  top = find(top_f);
else 
  wall = find(wall_f);
  top = find(top_f & ! wall_f);
endif

asave("cubcav.nod.tmp",x3);
asave("cubcav.con.tmp",i3);
if use_prismatic
  ico_prism = [i3(:,[1 2 3 5 6 7]);
	       i3(:,[3 4 1 7 8 5])];
  asave("cubcav.con-prism.tmp",ico_prism);
elseif use_tetra
  system("../../tools/hexasplit.bin -i cubcav.con.tmp -o cubcav.con-tetra.tmp");
endif


fid = fopen("cubcav.fixa.tmp","w");
fid != -1 || error("Couldn't open cubcav.fixa.tmp");

for j=wall'
  if strcmp(CASE,"laplace")
    fprintf(fid,"%d 1 0.\n",j);
  else
    fprintf(fid,"%d 1 0.\n%d 2 0.\n%d 3 0.\n",j,j,j);
  endif
endfor

for j=top'
  if strcmp(CASE,"laplace")
    fprintf(fid,"%d 1 1.\n",j);
  else
    fprintf(fid,"%d 1 1.\n%d 2 0.\n%d 3 0.\n",j,j,j);
  endif
endfor
fclose(fid);

if strcmp(case_in,"srfgath")
  fid = fopen("cubcav.fixa-srfgath.tmp","w");
  for j=create_set([top' wall']);
    for k=1:4
      fprintf(fid,"%d 1 %.20f\n%d 2 %.20f\n%d 3 %.20f\n%d 4 %.20f\n",
	      j,x3(j,1),j,x3(j,2),j,x3(j,3),j,x3(j,1));
    endfor
  endfor
  fclose(fid);
endif

