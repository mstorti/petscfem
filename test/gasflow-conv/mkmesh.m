###key mkmesh.m
### $Id: $

source("data.m.tmp");

w = zhomo([0,L/n,0,L],2,n+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2,1]);

asave("gfnwttest.nod.tmp",xnod);
asave("gfnwttest.con.tmp",icone);

pfperi("gfnwttest.peri.tmp",n+1+(1:n+1)',(1:n+1)',1:4);
pfperi("gfnwttest.peri-x.tmp",n+1,1,1:4);

u1 = [rhoref,uref,pref];
u2 = compshock(u1);

u1 = u1(:,[1,2,3,3]);
u1(3)=0;
u2 = u2(:,[1,2,3,3]);
u2(3)=0;

pffixa("gfnwttest.fixa.tmp",1,1:4,u1);
# pffixa("gfnwttest.fixa.tmp",n+1,1:4,u2,"a");
# pffixa("gfnwttest.fixa.tmp",n+1,4,u2(4),"a");

nnod = rows(xnod);
uini1 = u1;
uini1 = uini1(ones(nnod,1),:);

uini2 = u2;
uini2 = uini2(ones(nnod,1),:);
## xi = (xnod(:,1)>L/2);
## This is smoother
xi = (1+tanh((xnod(:,1)-L/4)/0.1))/2;
xi = xi - (1+tanh((xnod(:,1)-3*L/4)/0.1))/2;
## This is periodic in x
xi = xnod(:,1)/L;
xi = 0.5*(1-tanh(7*cos(2*pi*xi)));

uini = uini1+leftscal(xi,(uini2-uini1));
uini = uini1;

uini = [uini;zeros(1,4)];
asave("gfnwttest.ini.tmp",uini);

xnod = [xnod;L,0];
asave("gfnwttest.nod.tmp",xnod);
