u=[];
v=[];
w=[];
p=[];
rho=[];
k=0;
N = Nr+1;
while 1
  file = ["strip.state_" int2str(k) ".tmp"];
  [info, err, msg] = lstat (file);
  if err; break; endif
  eval(["U=aload(\"" file "\");"]);
  rho = [rho U(1:N,1)];
  u = [u U(1:N,2)];
  v = [v U(1:N,3)];
  w = [w U(1:N,4)];
  p = [p U(1:N,5)];
  k = k+1;
endwhile
printf("%d states loaded\n",k);

ga = 1.4;
g1 = ga-1;
Rgas = 287;
T = p./rho/Rgas;
ene = p./rho/g1;
v2 = u.^2+v.^2+w.^2;
Ene = ene + 0.5*v2;
rEne = rho.*Ene;

rene_int = sum(leftscal(diff(x).*(2*pi*xcent(x)),xcent(rho.*ene)))';
rv2_int = sum(leftscal(diff(x).*(2*pi*xcent(x)),xcent(0.5*rho.*v2)))';
rEne_int = sum(leftscal(diff(x).*(2*pi*xcent(x)),xcent(rEne)))';
