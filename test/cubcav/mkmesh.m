source("data.m.tmp");
#  N=5;
#  hratio = 0;
#  tol=1e-3/N;
#  leaky = 0;

w=zhomo([0 1 0 1],N+1,N+1,[1 hratio 1 1 hratio 1]);
[x2,i2] = pfcm2fem(w);
x2 = x2(:,[2 1]);

[x3,i3] = extrude(x2,i2,N,1/N);

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

fid = fopen("cubcav.fixa.tmp","w");
fid != -1 || error("Couldn't open cubcav.fixa.tmp");

for j=wall'
  fprintf(fid,"%d 1 0.\n%d 2 0.\n%d 3 0.\n",j,j,j);
endfor

for j=top'
  fprintf(fid,"%d 1 1.\n%d 2 0.\n%d 3 0.\n",j,j,j);
endfor
fclose(fid);
