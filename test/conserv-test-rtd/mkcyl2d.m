## __INSERT_LICENSE__
## $Id: mkductadv.m,v 1.18 2004/09/06 00:01:45 mstorti Exp $

source("./data.m.tmp");

Vmax = 1.5*Vmean;

## build a square
w=zhomo([0,Lz,-radius,radius],nz+1,n+1,[1 0 1 1 rratio 1]);
[xnod,icone]=pfcm2fem(w);
icone = icone(:,[1,4,3,2]);
nnod=rows(xnod);

asave("ductadv.con.tmp",icone);
asave("ductadv.nod.tmp",xnod);

inlet = (1:n+1)';
outlet = (n+1)*nz+(1:n+1)';
wall = (0:nz)'*(n+1);
wall = [1+wall;
        n+1+wall];

## bcconv connectivity
bcc = [outlet(2:n+1),outlet(1:n)];
asave("ductadv.bcconv-outlet-con.tmp",bcc);

tmp = complement(inlet,wall);
pffixa("ductadv.fixa-wall.tmp",tmp,(1:3));

y = xnod(inlet,2);
uinlet = [(1-(y/radius).^2)*Vmax,0*y];
pffixa("ductadv.fixa-inlet.tmp",inlet,1:2,uinlet);
  
tmp = complement(wall,outlet);
pffixa("ductadv.fixa-u-outlet.tmp",outlet,2);

pffixa("ductadv.fixa-inlet.tmp",inlet,1,1.0);

pfperi("ductadv.peri.tmp",inlet,outlet,1:3);

uini = [Vmean,0,0];
nnod = rows(xnod);
uini = uini(ones(nnod,1),:);
asave("ductadv.ini.tmp",uini);

## For the conservation test
if conserv_test
  ## u= a Dirac's delta in the center
  n1 = round(nz/2);
  nnod = rows(xnod);
  u = zeros(nnod,1);
  u(n1*(n+1)+(1:n+1)') = 1;
  asave("ductadv.ini.tmp",u);
  
  pfperi("ductadv.peri.tmp",inlet,outlet,1);
endif

ux = (1-(xnod(:,2)/radius).^2)*Vmax;
asave("ductadv.nod-vel.tmp",[xnod,0*ux,ux]);
