##__INSERT_LICENSE__
## $Id: mkcondwall.m,v 1.2 2005/03/28 16:42:59 mstorti Exp $

source("data.m.tmp");

w1 = zhomo([0 Lx1 0 Ly],Nx1+1,Ny+1,[hratio1 0 1 1 0 1]);
[xnod1,icone1] = pfcm2fem(w1);
## xnod1 = xnod1(:,[2 1]);
nnod1 = rows(xnod1);

w2 = zhomo([Lx1 Lx1+Lx2 0 Ly],Nx2+1,Ny+1,[1 0 hratio2 1 0 1]);
[xnod2,icone2] = pfcm2fem(w2);
## xnod2 = xnod2(:,[2 1]);
nnod2 = rows(xnod2);

xnod = [xnod1;xnod2];
icone = [icone1;
	 icone2+nnod1];
icone = icone(:,[1 4 3 2]);

## gplfem(xnod,icone);

left1 = (1:Ny+1)';
right1 = Nx1*(Ny+1)+(1:Ny+1)';

left2 = nnod1+(1:Ny+1)';
right2 = nnod1+Nx2*(Ny+1)+(1:Ny+1)';

base1 = (0:Nx1)'*(Ny+1)+1;
top1 = base1+Ny;

base2 = nnod1+(0:Nx2)'*(Ny+1)+1;
top2 = base2+Ny;

## Periodic in y direction
ndof = 3;
pfperi("condwall.peri1.tmp",base1,top1,(1:ndof));
pfperi("condwall.peri2.tmp",base2,top2,(1:ndof));

## Pressure at inlet
tmp = complement(top1,left1);
pffixa("condwall.fixa-in.tmp",tmp,ndof,DP);

## No flux at the internal wall (left side)
tmp = complement(top1,right1);
pffixa("condwall.fixa-wall1.tmp",tmp,1:2);

## No flux at the internal wall (rightt side)
tmp = complement(top2,left2);
pffixa("condwall.fixa-wall2.tmp",tmp,1:2);

## Pressure at outlet
tmp = complement(top2,right2);
pffixa("condwall.fixa-out.tmp",tmp,ndof,0);

## Two layers of fictitious nodes at the internal wall
fic1_base = nnod1+nnod2;
fic1 = fic1_base+(1:Ny+1);

xnod = [xnod;
	xnod(right1,:)];

## Two layers of fictitious nodes at the internal wall
fic1_base = nnod1+nnod2;
fic1 = fic1_base+(1:Ny+1)';

xnod = [xnod;
	xnod(left2,:)];

fic2_base = nnod1+nnod2+(Ny+1);
fic2 = fic2_base+(1:Ny+1)';

## Connectivities for the `cond_wall' elemeset
icowall = [right1,left2,fic1,fic2];
## icowall(Ny+1,:)=[]; ## Last node goes by periodic b.c.'s

asave("condwall.nod.tmp",xnod);
asave("condwall.con.tmp",icone);
asave("condwall.condwall-con.tmp",icowall);

## Connectivities for the `cond_wall' elemeset
pfperi("condwall.wall-peri.tmp",right1,left2,(1:ndof)');
