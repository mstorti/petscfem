N=128;
volume_relax_coef = 0.1;
nfilt = 20;

if 0
  u = zeros(N,1);
  u(round(N/2))=1;
elseif 0
  u = -cos(6*2*pi*(0:N-1)'/N);
elseif 1
  u = (1:N)';
  u = (u>=(N/2));
endif
u0 = u;

indx1= (0:N-1)';
indx0 = modulo(indx1-1,N)+1;
indx2 = modulo(indx1+1,N)+1;
indx1 = indx1+1;

for k=1:nfilt
  w = u + volume_relax_coef * (u(indx2) - 2*u + u(indx0));
  u = w;
endfor
