source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

xnod = [xnod eta0+(etaL-eta0)/Lx*xnod(:,1)];

asave("aqui.nod.tmp",xnod);
asave("aqui.con.tmp",icone);

x=xnod(:,1);
y=xnod(:,2);

tol=1e-5;

fid = fopen("aqui.fixa.tmp","w");
x0 = (1:Ny+1)';
for k=x0'
  fprintf(fid,"%d 1 0.\n",k);
endfor

xL = Nx*(Ny+1)+(1:Ny+1)';
for k=xL'
  fprintf(fid,"%d 1 1.\n",k);
endfor
fclose(fid);
