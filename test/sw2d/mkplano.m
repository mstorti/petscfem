source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 0 1 1 0 1]);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);
nnod=rows(xnod);

zb=zeros(nnod,1);
xnod=[xnod(:,1:2) zb];
asave("plano.nod.tmp",xnod);
asave("plano.con.tmp",icone);

for j=1:nnod
  h_ini(j) = ((h_out-h_in)/Lx).*xnod(j,1)+h_in;
endfor
if (turbulent)
  ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_ini*ones(nnod,1),\
       kappa*ones(nnod,1),epsilon*ones(nnod,1)];
else
##  ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_ini*ones(nnod,1)];
  ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_ini'];
endif
asave("plano.ini.tmp",ini);

tol=1.e-4;
in=(1:Ny+1)';
out=((Nx)*(Ny+1)+1:(Nx+1)*(Ny+1))';
## l=(2*(Ny+1):Ny+1:(Nx)*(Ny+1))';
## r=((Ny+1)+1:Ny+1:(Nx)*(Ny+1))';
l=((Ny+1):Ny+1:(Nx+1)*(Ny+1))';
r=(1:Ny+1:(Nx+1)*(Ny+1))';

##pffixa("plano.in-fixa.tmp",in,[1 2],[u_in,v_in]);
##pffixa("plano.out-fixa.tmp",out,[3],[h_out]);
pffixa("plano.in-fixa.tmp",in,[3],[h_in]);
pffixa("plano.out-fixa.tmp",out,[1,2],[u_out,v_out]);

if (turbulent)
  pffixa("plano.ke-fixa.tmp",(1:nnod)',[4 5],[kappa,epsilon]);
endif
pfconstr("plano.constr-slip-l.tmp",[l,l],1:2,[0. 1.]);
pfconstr("plano.constr-slip-r.tmp",[r,r],1:2,[0. -1.]);
pfperi("plano.peri-lr.tmp",l,r,[1:3]);
pffixa("plano.v-fixa-all.tmp",l,[2],0);
pffixa("plano.h-fixa-in.tmp",in,[3],h_in);
