source("data.m.tmp");

hav=L/N;
w=zhomo([0 L 0 Ly],N+1,Ny+1,[1 0 hratio 1 0 hratioy]);
[xnod,icone]=pfcm2fem(w);
xnod = [xnod;
        xnod(1:Ny+1,:)];

icone=icone(:,[1 4 3 2]);
nonlr=[(1:Ny+1)',(N+1)*(Ny+1)+(1:Ny+1)'];
wallke = [(1:Ny)' (2:Ny+1)'];

fid = fopen("wallke.fixa_lag.tmp","w");
for k=(N+1)*(Ny+1)+(1:Ny+1)
  for j=3:5
    fprintf(fid,"%d %d %f\n",k,j,0.);
  endfor
endfor
fclose(fid);

uin=[0.66 0. 0. 2e-4  1e-05];  
fid = fopen("wallke.fixa_in.tmp","w");
for k=1:(Ny+1):(N+1)*(Ny+1)
  for j=1:5
    fprintf(fid,"%d %d %f\n",k,j,uin(j));
  endfor
endfor
fclose(fid);

inlet = (1:(Ny+1):(N+1)*(Ny+1))';
outlet = (Ny+1:(Ny+1):(N+1)*(Ny+1))';

fid = fopen("wallke.peri.tmp","w");
for k=1:N+2
  in = (k-1)*(Ny+1)+1;
  out = k*(Ny+1);
  for j=1:5
    fprintf(fid,"%f %d %d    %f %d %d \n",-1.,out,j,+1.,in,j);
  endfor
endfor
fclose(fid);

uini=[0. 0.666 0. 0.1 0.05];
uini=[uini(ones((N+1)*(Ny+1),1),:);
      zeros(Ny+1,5)];

boux0=[(1:Ny+1)' ones(Ny+1,1) zeros(Ny+1,1)];

asave("wallke.nod.tmp",xnod);
asave("wallke.con.tmp",icone);
asave("wallke.ini.tmp",uini);
asave("wallke.nonlr.tmp",nonlr);
asave("wallke.wallke.tmp",wallke);
asave("wallke.u0.tmp",boux0);
