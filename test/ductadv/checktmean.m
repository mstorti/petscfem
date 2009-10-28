###key proc.m
### $Id: $

source("data.m.tmp");

R1 = 0.*radius;

if !exist("xnod");
  xnod = aload("ductadv.nod-vel.tmp");
  ux = xnod(:,3);
  xnod = xnod(:,1:2);
  icone = aload("ductadv.con.tmp");
endif

nnod = rows(xnod);
rem(nnod,nz+1)==0 || error("bad nnod");
nlay = nnod/(nz+1);
outlet = (n+1)*nz+(1:n+1)';

kk = 1:2;
if ndim==2; kk=2; endif
r = l2(xnod(1:nlay,kk));
d = abs(r-R1);
[bid,ix] = min(d);
clear bid
indx = ix + (0:nz)'*nlay;

kinc = 1;
kend = -1;
dire = "./STEPS";
kk = 3;
if ndim==2; kk=1; endif

if !exist("uh");
  su = [];
  uh = [];
  cout = [];
  k = 0;
else 
  k = columns(uh);
endif

nread = 0;
while 1
  file = sprintf("%s/ductadv.state-%d.tmp.gz",dire,k);
  if !exist(file); break; endif
  # printf("processing file %s\n",file);
  u = aload(file);
  nread++;
  uh = [uh,u(indx)];
  su = [su;pfint(xnod,icone,u)];
  icout = pfint(xnod(:,2),outlet,u);
  icoutu = pfint(xnod(:,2),outlet,u.*ux);
  cout = [cout;icout,icoutu];
  k += kinc;
  if kend>=0 && k>=kend; break; endif
endwhile
printf("read %d states\n",nread);

#tmean = sum(max(cout(:,2))-cout(:,2))*kinc*Dt/(Vmean*2*radius);
Qout = pfint(xnod(:,2),outlet,ux);
tmean = sum(1-cout(:,2)/Qout)*kinc*Dt;
tmeanv = cumsum(1-cout(:,2)/Qout)*kinc*Dt;

tmeananal = Lz/Vmean;
printf("tmean computed %f, analytic %f\n",tmean,tmeananal);

tol = 0.15;
erro = abs(tmean-tmeananal);
relerro = erro/tmeananal;
printf("OK ? %d (abs. error %f, rel error %f, rel.err. tol %f)\n", \
       relerro<tol,erro,relerro,tol);
