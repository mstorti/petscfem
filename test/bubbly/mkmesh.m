source("data.m.tmp");

rem(N,2)==0 || warning("N should be even");

w=zhomo([0 1 0 1],N+1,N+1,[1 hratio 1 1 hratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

asave("bubbly.nod.tmp",xnod);
asave("bubbly.con.tmp",icone);

x=xnod(:,1);
y=xnod(:,2);

tol=1e-5;
lid=find(abs(y-1)<tol)';
nlid=length(lid);
b=create_set([find(abs(x)<tol);
              find(abs(x-1)<tol);
              find(abs(y)<tol)]);
b=complement(lid,b);
nb=length(b);

if g_body 
  utop = 0;
else
  utop = 1;
endif

ul=3;
vl=4;
fixa=[lid' ones(nlid,2)*diag([ul utop]);
      lid' ones(nlid,2)*diag([vl 0]);
      b' ones(nb,2)*diag([ul 0]);
      b' ones(nb,2)*diag([vl 0])];

asave("bubbly.fixa.tmp",fixa);
fid = fopen("bubbly.fixa_all.tmp","w");
for j=1:(N+1)^2;
  fprintf(fid,"%d %d %f\n",j,1,alpha_l);
  fprintf(fid,"%d %d %f\n",j,5,u_g);
  fprintf(fid,"%d %d %f\n",j,6,v_g);
  fprintf(fid,"%d %d %f\n",j,7,k);
  fprintf(fid,"%d %d %f\n",j,8,eps);
endfor
fclose(fid);

## fix all quantities not liquid velocity and pressure


## y coordinates of nodes on the centerline
ny = (N/2)*(N+1)+(1:N+1)';
yh = xnod(ny,2);

save bubbly.ny.tmp  ny yh
