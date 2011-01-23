## $Id: mkmesh.m,v 1.3 2006/03/12 12:07:37 mstorti Exp $
source("data.m.tmp");

# Ntrail = 20;
# phi = (0:Ntrail)'/Ntrail*pi;
# xtra = [cos(phi),2.5*ones(size(phi)),sin(phi)];
# icotra = [(1:Ntrail)',(1:Ntrail)'+1];
# mkvtkfile("trail.vtk","hoho",xtra,icotra,xtra(:,1),"xx");

w = zhomo([0 Ly 0.5*[Rratio*Lx Lx]],Ny+1,Nx+1);
[x2,ico2] = pfcm2fem(w);
x2 = x2(:,[2,1]);

[x0,icone] = extrude(x2,ico2,Nphi,1/Nphi);
y = x0(:,2);
r = x0(:,1);
Lrfacy = 0.4;
# rfacy0 = 1.5;
rfacy0 = 1.;
rfacy = 1+(rfacy0-1)./(1+y/Lrfacy);
r = r.*rfacy;
phi = 2*pi*x0(:,3);
x0 = [r.*cos(phi),y,r.*sin(phi)];
eta = y/Ly;
# eta -= 0.75*eta.*(1-eta);
y = eta*Ly;
x0(:,2) = y;

nlay = (Nx+1)*(Ny+1);
nnod = rows(x0);
nnod = Nphi*(nlay+1) || error("bad nnod");
select = (1:nlay)';
select = [select;select+Nphi*nlay];
[x0,icone,iren]=pasten(x0,icone,select);
[x0,icone] = pfextmesh(x0,icone);
nnod = rows(x0);

asave("elastld.nod.tmp",x0);
asave("elastld.ini.tmp",zeros(nnod,3));
asave("elastld.con.tmp",icone);

vel = 1;
xold = zeros(size(x0));
xold(:,1) = -vel*Dt;
asave("elastld.iniold.tmp",xold);

tol = 1e-5;
bot = find(x0(:,2)<tol);
pffixa("elastld.fixa-bot.tmp",bot,1:ndim,1);

mkvtkfile("beam.vtk","haha",x0,icone,x0(:,1),"x");

[bid,indx] = max(x0(:,2));
asave("elastld-nodes.tmp",indx);

## Save data for computation of the trail
nodes = (0:Nphi-1)'*nlay+(Nx+1)*Ny+1;
fid = fopen("elastld.trail.tmp","w");
for k=1:length(nodes)
  fprintf(fid,"%d %f\n",nodes(k),1.0);
endfor
fclose(fid);
