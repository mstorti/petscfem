## $Id: mkplate.m,v 1.6 2005/01/31 02:20:52 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);

uini = Machin*cref;
poutlet = pref;

w = zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 0 1 1 0 yratio]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

indx = find(xnod(:,1)<Lslip);
x1 = [Lslip;xnod(indx,1)];
x1 = onedstr([10 0 1],x1);
x1(1) = [];
xnod(indx,1) = x1;

indx = find(xnod(:,1)>=Lslip);
x1 = [Lslip;xnod(indx,1)];
x1 = onedstr([1 0 10],x1);
x1(1) = [];
xnod(indx,1) = x1;

tol=1e-5;
inlet = find(abs(xnod(:,1))<tol);
outlet = find(abs(xnod(:,1)-Lx)<tol);

wall = find(abs(xnod(:,2))<tol);
slip = find(xnod(wall,1)<Lslip \
	    || xnod(wall,1)>Lslip+Lplate);
slip = wall(slip);
wall = complement(slip,wall);
outer = find(abs(xnod(:,2)-Ly)<tol);

## rho,u,v at inlet
pffixa("gfplate.fixa-in.tmp", \
       inlet,1:3,[rhoref uini 0.0])

## v=0 at slip
tmp = complement(inlet,slip);
pffixa("gfplate.fixa-slip.tmp",slip,3);

## u=v=0 at wall
pffixa("gfplate.fixa-wall.tmp",wall,[2 3]);

## v=0 at outer
tmp = complement(inlet,outer);
pffixa("gfplate.fixa-outer.tmp",tmp,3);

## p fixed at outlet
## outlet = complement(wall,outlet)';
pffixa("gfplate.fixa-outlet.tmp", \
       outlet,4,pref);

nnod = size(xnod,1);

## Absorbing b.c.'s at outlet
ficabso = nnod+(1:length(outlet))';
nnod2 = nnod+length(outlet);

asave("gfplate.abso-con.tmp",[outlet,ficabso,ficabso]);
xnod = [xnod;
	xnod(outlet,:)];

## Fictitious nodes for fixed
## wall temperature b.c.
ficwall = nnod2+(1:Nx+1)';
nnod2 = nnod2+Nx+1;

xnod = [xnod;
	zeros(Nx+1,2)];
fid = fopen("gfplate.twall.tmp","w");
for k=1:Nx+1
  node = (k-1)*(Ny+1)+1;
  ficnode = ficwall(k);
  fprintf(fid,"%d %d\n",node,ficnode);
  xnod(ficnode,:) = xnod(node,:);
endfor
fclose(fid);

pffixa("gfplate.fixa-twall.tmp", \
       ficwall,2:4,[Tref,0,0]);

## This is used if not Twall imposed
pffixa("gfplate.fixa-lagmul-tw.tmp", \
       ficwall,1);

asave("gfplate.nod.tmp",xnod);
asave("gfplate.dx-nod.tmp",xnod(1:nnod,:));
asave("gfplate.con.tmp",icone);

Uini = [rhoref,uini,0,pref];
Uini = Uini(ones(nnod2,1),:);
Uini(ficwall,:) = 0;
Uini(ficabso,:) = 0;
asave("gfplate.ini.tmp",Uini);
