## $Id: mkgfshock.m,v 1.8 2005/02/27 12:59:52 mstorti Exp $
source("data.m.tmp");

## Tuyere data
plt_file = "~/PETSC/COMP-CORNER/Z23_contour.plt";
z = aload(plt_file);
xtuy = z(:,1);
Atuy = z(:,2).^2;
Atuy_char = sqrt(max(Atuy)*min(Atuy));
x0 = min(z(:,1));
x1 = max(z(:,1));
Ltuy = x1-x0;

## Make width of strip proportional to `tuyere' Area
## Scale so that geometrical mean will be O(h)
ytuy = Atuy/Atuy_char*Ltuy/Nx;

## Make 1D mesh
w = zhomo([0 Lx/Nx 0 Lx ],2,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);
nnod = size(xnod,1);
x = xnod(1:nnod,1);
x = x0+(x/Lx)*Ltuy;

## Map sides to tuyere area and compute derivatives
## in order to impose slip condition
xx = x(1:Nx+1,1);
yy = spline(xtuy,ytuy,xx);
dx = 1e-4*Ltuy;
yyp = spline(xtuy,ytuy,xx+dx);
yym = spline(xtuy,ytuy,xx-dx);
dydx = (yyp-yym)/(2*dx);

xnod((1:Nx+1)',2) = -yy;
xnod(Nx+1+(1:Nx+1)',2) = +yy;

asave("gfshock.nod.tmp",xnod);
asave("gfshock.con.tmp",icone);

## rho,u,v at inlet
inlet = [1];
pffixa("gfshock.fixa-in.tmp", \
       inlet,[1 3 4],[rhoin0 0 pin0]);

## p at outlet
outlet = [Nx+1];
pffixa("gfshock.fixa-outlet.tmp", \
       outlet,4,pout0)

## Symmetry in y for rho, u, p
pfperi("gfshock.peri.tmp", \
       Nx+1+(1:Nx+1)',(1:Nx+1)',[1 2 4]);

## Antisymmetry for v
pfconstr("gfshock.v-peri.tmp", \
	 [Nx+1+(1:Nx+1)',(1:Nx+1)'],3,[1 1]);

## Slip on the low border
pfconstr("gfshock.slip.tmp", \
	 [(1:Nx+1)',(1:Nx+1)'],[2 3], \
	 [dydx,ones(Nx+1,1)]);

## Fixa on reference nodes
Uini = [rhoout0,0,0,pout0];
Uini = Uini(ones(nnod,1),:);

asave("gfshock.ini.tmp",Uini);
asave("gfshock.some-nodes.tmp",(1:Nx+1)');
