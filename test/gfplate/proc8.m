## Check conservativity
source("data.m.tmp");
## resdir = "/cdrom/gfplate-results/gfplate-2005-FEB-26e";
resdir = ".";

axi = 1;
if 1
  u0 = aload([resdir "/gfshock3d.some-rslt.tmp"]);
  X = aload([resdir "/gfshock3d.nod.tmp"]);
endif
u = u0;

nodes_some = create_set(u(:,1))';
nsome = length(nodes_some);
all(nodes_some==u(1:nsome,1)) || error("malformed `some' file");
u(:,1) = [];

if !axi
  x = X(nodes_some,1);
  A = -2*X(nodes_some,2);
  rhoindx = 1;
  uindx = 2;
  pindx = 4;
else
  x = X(nodes_some,3);
  A = pi*X(nsome+(1:nsome),1).^2;
  rhoindx = 1;
  uindx = 4;
  pindx = 5;
endif

nt = rows(u);
rem(nt,nsome) ==0 || error("not correct number of rows in some file");
nt = nt/nsome;
rho = reshape(u(:,rhoindx),nsome,nt);
uu = reshape(u(:,uindx),nsome,nt);
p = reshape(u(:,pindx),nsome,nt);
c2 = ga*p./rho;

G = leftscal(A,rho.*uu);
Gmom = leftscal(A,rho.*uu.*uu+p);
Gene = leftscal(A,rho.*uu.*(0.5*uu.^2+c2/(ga-1)));

title("Flow mass conservation");
plot(x,G(:,1:10:nt));
pause

title("Flow mass conservation - last state");
plot(x,G(:,nt));
pause

title("Entalpy conservation");
plot(x,Gene(:,1:10:nt));
pause

title("Entalpy conservation - last state");
plot(x,Gene(:,nt));
