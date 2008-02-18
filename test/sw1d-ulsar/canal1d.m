#$Id: canal1d.m,v 1.4 2004/07/15 18:33:23 rodrigop Exp $
#crea malla y condiciones de borde e iniciales  de sw-1d
#condiciones absorbentes
## Number os elements in x, acordarse de poner nodos dummy
clear;
source("data.m.tmp");
slope=0.0; % 10% de la long del rio
x=Lx*onedstr([.1 .1 .1],Nx+1);# n de nodos deberia ser impar para que el bump quede en el centro
nnod=length(x);
y=zeros(nnod,1);
z=(-x*sqrt(2)*4e-5)*slope;
xnod=[x y z];
## add fictitious nodes
xnod=[xnod;xnod(1,:);xnod(nnod,:)];
icone1=[1:nnod-1]';
icone2=[2:nnod]';
icone=[icone1 icone2];
asave("canal1d.nod.tmp",xnod);
asave("canal1d.con.tmp",icone);
##parabolic bump paper mario for initial state
xc=xnod(round(Nx/2),1);
for i=1:nnod,
  h_iniv(i)=1.+A*((1./(sqrt(2.*pi))).*exp(-0.5*((xnod(i,1)-mu)/sig).^2));
endfor
h_iniv(nnod+1)=0.;h_iniv(nnod+2)=0.;
h_iniv=h_iniv';
u0v=u_ini*ones(nnod,1);
u0v(nnod+1)=0.;u0v(nnod+2)=0.;

asave("canal1d.ini.tmp",[u0v h_iniv]);
##pffixa("canal1d.in-fixa.tmp",[1],[1 2],[u_in,h_in]);
##pffixa("canal1d.out-fixa.tmp",[nnod],[1 2],[u_out,h_out]);
pffixa("canal1d.in-fixa.tmp",[1],[1],[u_in]);
pffixa("canal1d.out-fixa.tmp",[nnod],[2],[h_out]);

if (uref==0)
  asave("canal1d.abso-out.tmp",[1 nnod+1 -1.0; \
				nnod nnod+2 1.0]);
else
  asave("canal1d.abso-out.tmp",[1 nnod+1 nnod+1 -1.0;\
				nnod nnod+2 nnod+2  1.0]);
endif
