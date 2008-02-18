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
A=2.; 
for i=1:nnod,
h_iniv(i)=1.+A.*((1/(sqrt(2.*pi)))*exp(-0.5.*((xnod(i,1)-mu)/sig).^2));
end
h_iniv(1)=1.;h_iniv(nnod)=1.;

if (turbulent)
  ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_iniv',\
       kappa*ones(nnod,1),epsilon*ones(nnod,1)];
else
  ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_iniv*ones(nnod,1)];
endif
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

# pffixa("plano.in-fixa.tmp",in,[1 3],[u_in,h_in]);
# pffixa("plano.out-fixa.tmp",out,[3],[h_out]);
pffixa("plano.in-fixa.tmp",in,[1 2],[u_in,v_in]);
pffixa("plano.out-fixa.tmp",out,[1 2],[u_out,v_out]);
if (turbulent)
  pffixa("plano.ke-fixa.tmp",(1:nnod)',[4 5],[kappa,epsilon]);
endif

pfperi("plano.constr-wall-peri.tmp",l,r,1:3);

if (ulsar)
  asave("plano.abso-out.tmp",[2 nnod+1 -1.0 0.0; nnod nnod+2 1.0 0.0]);
else
  asave("plano.abso-out.tmp",[2 nnod+1 nnod+1 -1.0 0.0; nnod nnod+2 \
			      nnod+2 1.0 0.0]);
endif

## pfconstr("plano.constr-slip-l.tmp",[l,l],1:2,[0. 1.]);
## pfconstr("plano.constr-slip-r.tmp",[r,r],1:2,[0. -1.]);
