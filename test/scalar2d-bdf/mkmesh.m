## $Id: mkmesh.m,v 1.8 2007/01/28 00:29:19 mstorti Exp $
source("data.m.tmp");

w = zhomo([0 Lx 0 Ly],N+1,N+1);
[xnod,icone]=pfcm2fem(w);
icone = icone(:,[1 4 3 2]);
icone = [icone(:,[1,2,3]);
         icone(:,[3,4,1])];

asave("gascont.nod.tmp",[xnod,xnod,xnod]);

xele = pfnd2ele(xnod,icone,xnod);
nelem = rows(icone);
vele = [0,0];
vele = vele(ones(nelem,1),:);

if strcmp(adv_case,"gaussian");
  xele -= 0.5;
  vele = [-xele(:,2),xele(:,1)];
elseif strcmp(adv_case,"gaussian_diag");
  vele = ones(nelem,2);
endif
asavecon("gascont.con.tmp",icone,vele);

nnod = rows(xnod);
mkvtkfile("gascont.vtk","gascont.vtk",xnod(:,1:2),ones(nnod,1),"u");

uini = ones(nnod,1);
if strcmp(adv_case,"gaussian");
  xx = xnod(:,1);
  yy = xnod(:,2);
  r = l2([xx-0.5*Ly,yy-0.25*Ly]);
  sigma = 0.1;
  uini = 1+exp(-(r/sigma).^2);
elseif strcmp(adv_case,"gaussian_diag");
  xx = xnod(:,1);
  yy = xnod(:,2);
  r = l2([xx-0.25*Ly,yy-0.25*Ly]);
  sigma = 0.1;
  uini = 1+exp(-(r/sigma).^2);
endif

asave("gascont.ini.tmp",uini);

tol = 1e-5;
bot = find(abs(xnod(:,2))<tol);
top = find(abs(xnod(:,2)-1)<tol);
left = find(abs(xnod(:,1))<tol);
right = find(abs(xnod(:,1)-1)<tol);
tmp = unique([bot;top;left;right]);
pffixa("gascont.fixa.tmp",tmp,1,1);

pfperi("gascont.periy.tmp",top,bot,1)
pfperi("gascont.perix.tmp",right,left,1)

# asave("gascont.some-nodes.tmp",(1:Ny+1)');
