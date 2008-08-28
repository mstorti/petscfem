source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 0 1 1 0 1]);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);
nnod=rows(xnod);

zb=zeros(nnod,1);
xnod=[xnod(:,1:2) zb];
## add fictitious node 
xnod=[xnod;xnod(2,:);xnod(nnod,:)];
asave("plano.nod.tmp",xnod);
asave("plano.con.tmp",icone);

xc=xnod(round(Nx/2),1);
mu=Lx/2;
sig=0.04*Lx;
rc=.1*Lx;
A=.1;
for i=1:nnod,
h2_iniv(i)=1.+A.*((1/(sqrt(2.*pi)))*exp(-0.5.*((xnod(i,1)-mu)/sig).^2));
h1_iniv(i)=h1_ini;
end
h1_iniv(1)=h1_ini;h1_iniv(nnod)=h1_ini;
h2_iniv(1)=h2_ini;h2_iniv(nnod)=h2_ini;

ini=[u1_ini*ones(nnod,1),v1_ini*ones(nnod,1),h1_iniv',\
     u2_ini*ones(nnod,1),v2_ini*ones(nnod,1),h2_iniv'];

#ini=[u1_ini*ones(nnod,1),v1_ini*ones(nnod,1),h1_ini*ones(nnod,1),\
#     u2_ini*ones(nnod,1),v2_ini*ones(nnod,1),h2_ini*ones(nnod,1)];

## add ini to fictitious node
ini=[ini;zeros(1,ndof);zeros(1,ndof)];
asave("plano.ini.tmp",ini);

tol=1.e-4;
in=(1:Ny+1)';
out=((Nx)*(Ny+1)+1:(Nx+1)*(Ny+1))';
## l=(2*(Ny+1):Ny+1:(Nx)*(Ny+1))';
## r=((Ny+1)+1:Ny+1:(Nx)*(Ny+1))';
l=((Ny+1):Ny+1:(Nx+1)*(Ny+1))';
r=(1:Ny+1:(Nx+1)*(Ny+1))';

## pffixa("plano.out-fixa.tmp",out,[3 6],[h1_out,h2_out]);
pffixa("plano.in-fixa.tmp",in,[1 2 4 5],[u1_in,v1_in,u2_in,v2_in]);
pffixa("plano.out-fixa.tmp",out,[3 6],[h1_out,h2_out]);
pffixa("plano.v-fixa.tmp",l,[2 5],[0,0]);

pfperi("plano.in-out-peri.tmp",in,out,1:6);
pfperi("plano.constr-wall-peri.tmp",l,r,1:6);

if (ulsar)
  asave("plano.abso-out.tmp",[2 nnod+1 -1.0 0.0; nnod nnod+2 1.0 0.0]);
else
  asave("plano.abso-out.tmp",[2 nnod+1 nnod+1 -1.0 0.0; nnod nnod+2 \
			      nnod+2 1.0 0.0]);
endif
