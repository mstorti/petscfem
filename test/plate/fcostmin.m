omega=0.25;
xi=2;

a = [1/omega^2+xi/omega+1/2;
     -2/omega^2;
     1/omega^2-xi/omega+1/2];

b = [0.5+xi/omega 0 0.5-xi/omega]';

a0 = [a(2:3);b(3)];
fcost(a0);
fcostg(a0);

## [x,info] = fsolve('fcostg',a0);
eps = 1e-4;
a = a0;
ah = [];
gh = [];
fh = [];
iter = 0;
niter=250;
while 1
  g = fcostg(a);
  ah = [ah; a'];
  gh = [gh; g'];
  fh = [fh; fcost(a)];
  a = a - 1*g;
  ah = [ah; a'];
  iter = iter+1;
  if iter==niter; break; endif
endwhile

x = fcost(a,1);
a = x(1:3);
b = x(4:6);
sfilter
