## $Id: mkplate.m,v 1.3 2005/01/08 15:26:04 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);

uini = Machin*cref;
poutlet = pref;

Ly=2;
yratio = 10;
w = zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 0 1 1 0 yratio]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

tol=1e-5;
inlet = find(abs(xnod(:,1))<tol);
outlet = find(abs(xnod(:,1)-Lx)<tol);

wall = find(abs(xnod(:,2))<tol);
slip = find(xnod(wall,1)<Lplate);
slip = wall(slip);
wall = complement(slip,wall);
outer = find(abs(xnod(:,2)-Ly)<tol);

## rho,u,v at inlet
pffixa("gfplate.fixa-in.tmp",inlet,1:3,[rhoref uini 0.0])

## v=0 at slip
tmp = complement(inlet,slip);
pffixa("gfplate.fixa-slip.tmp",slip,3);

## u=v=0 at wall
pffixa("gfplate.fixa-wall.tmp",wall,[2 3]);

## v=0 at outer
tmp = complement(inlet,outer);
pffixa("gfplate.fixa-outer.tmp",tmp,3);

## p fixed at outlet
tmp = complement(wall,outlet);
pffixa("gfplate.fixa-outlet.tmp",tmp,[3,4],[0,pref]);

Uini = [rhoref,uini,0,pref];

nnod = size(xnod,1);
Uini = Uini(ones(nnod,1),:);
asave("gfplate.ini.tmp",Uini);

## Fictitious nodes
ficnodes = nnod+(1:Nx+1)';

xnod = [xnod;
	zeros(Nx+1,2)];
fid = fopen("gfplate.twall.tmp","w");
for k=1:Nx+1
  node = (k-1)*(Ny+1)+1;
  ficnode = nnod+k;
  fprintf(fid,"%d %d\n",node,ficnode);
  xnod(ficnode,:) = xnod(node,:);
endfor
fclose(fid);

pffixa("gfplate.fixa-twall.tmp",ficnodes,2:4,[Tref,0,0]);

asave("gfplate.nod.tmp",xnod);
asave("gfplate.con.tmp",icone);
