source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

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
  fprintf(fid,"%d %d %f\n",j,2,alpha_l);
#  fprintf(fid,"%d %d %f\n",j,3,0.);
#  fprintf(fid,"%d %d %f\n",j,4,0.);
  fprintf(fid,"%d %d %f\n",j,5,u_g);
  fprintf(fid,"%d %d %f\n",j,6,v_g);
  fprintf(fid,"%d %d %f\n",j,7,k);
  fprintf(fid,"%d %d %f\n",j,8,eps);
endfor
fclose(fid);

nnod = length(x);
## uini = [rand(nnod,1) alpha_l rand(nnod,2) 0. 0. 0.1 0.1];
uini = [0. alpha_l 0. 0. 0. 0. 0.1 0.1];
uini=uini(ones(rows(x),1),:);
asave("bubbly.ini.tmp",uini);

bottom = (1:(N+1):(N+1)^2)';
right = (N+1)*N+(1:N+1)';
top = ((N+1)^2:-(N+1):(N+1))';
left = ((N+1):-1:1)';

bcconv = [bottom(1:N) bottom(2:N+1);
	  right(1:N) right(2:N+1);
	  top(1:N) top(2:N+1);
	  left(1:N) left(2:N+1)];
asave("bubbly.bcconv.tmp",bcconv);

## fix all quantities not liquid velocity and pressure


## y coordinates of nodes on the centerline
ny = (N/2)*(N+1)+(1:N+1)';
yh = xnod(ny,2);

save bubbly.ny.tmp  ny yh
