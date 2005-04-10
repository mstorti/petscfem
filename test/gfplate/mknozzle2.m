## $Id: mknozzle2.m,v 1.6 2005/04/10 12:24:30 mstorti Exp $
source("data.m.tmp");

Nx = Nx1+Nx2;
ndof=4;

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);

uini = Machin*cref;
poutlet = pref;

w = zhomo([0 Lx1+Lx2 0 Ly],Nx1+Nx2+1,Ny+1,[1 0 1 1 0 yratio]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

x = xnod(:,1);
y = xnod(:,2);

## Make contraction in the [0,Lx1] range
indx = find(x<Lx1);
c = 1-(DLy/Ly)*(1-abs((x-Lx1/2)/(Lx1/2)));
y(indx) = y(indx).*c(indx);

## Make expansion
indx = find(x>=Lx1);

x(indx) = onedstr([1 0 rratio],x(indx));
r = (x(indx)-Lx1);
theta = pi/2*y(indx)/Ly;
dx = Lx1+r.*cos(theta)-x(indx);
dy = r.*sin(theta);

## Affect displacement by a blending function
## depending on the distance to the corner
R = 0.6*Ly;
r = l2([x(indx)-Lx1,y(indx)-Ly]);
cr = (1-r/R);
cr = cr.*(cr>0);
c = 1-(1-y(indx)/Ly).^0.5.*cr;

x(indx) = x(indx) + c.*dx;
y(indx) = y(indx) + c.*dy;

xnod = [x,y];

if 0
  gplfem(xnod,icone);
  gplot "malla.gpl"
  return
endif

nnodact = rows(xnod);

nodout = (Ny+1)*Nx+(1:Ny+1)';
nodfic = nnodact+(1:Ny+1)';
nnodfic = Ny+1;
xnod = [xnod;
	xnod(nodout,:)];

normalout = xnod(nodout,:);
normalout = leftscal(1./l2(normalout),normalout);
asave("gfnozzle2.nod.tmp",xnod);
asave("gfnozzle2.con.tmp",icone);

fid = fopen("gfnozzle2.abso-con.tmp","w");
for k=1:rows(nodout)
  fprintf(fid,"%d %d %d     %.12g %.12g\n",
	  nodout(k),nodfic(k),nodfic(k),
	  normalout(k,1), normalout(k,2));
endfor
fclose(fid);

# asave("gfnozzle2.abso-con.tmp", \
#       [nodout,nodfic,nodfic]);

inlet = (1:Ny+1)';
slip = 1+(Ny+1)*(0:Nx)';
outlet = Nx*(Ny+1)+inlet;
wall = (Ny+1)*(1:Nx+1);

## Tangents at the wall
nwall = length(wall);
tan = zeros(nwall,2);
## Tangents to elements
dtan = xnod(wall(2:nwall),:)-xnod(wall(1:nwall-1),:);
dtan = leftscal(1./l2(dtan),dtan);

tan(2:nwall-1,:) = 0.5*(dtan(1:nwall-2,:) \
			+dtan(2:nwall-1,:));
tan(1,:) = dtan(1,:);
tan(nwall,:) = dtan(nwall-1,:);
tan = leftscal(1./l2(tan),tan);

## Normals
nor = [-tan(:,2),tan(:,1)];

xwall = xnod(wall,:);
xwall2 = xwall+0.1*nor;

## Slip at the wall
tmp = complement(inlet,wall)';
pfconstr("gfnozzle2.constr-wall.tmp", \
	 [tmp,tmp],2:3,nor);

## [rho v p] at inlet
pffixa("gfnozzle2.fixa-in.tmp",inlet, \
       [1 3 4],[rhoref,0,pratio*pref]);

## `p' at outlet
pffixa("gfnozzle2.fixa-out.tmp",outlet,4,pref);

## `v=0' at axis
pffixa("gfnozzle2.fixa-slip.tmp",slip,3);

nnod = rows(xnod);
uini = [rhoref,0,0,pref];
uini = [uini(ones(nnodact,1),:);
	zeros(nnodfic,ndof)];
asave("gfnozzle2.ini.tmp",uini);
