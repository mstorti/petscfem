#$Id: proc3.m,v 1.1 2003/06/06 17:36:02 mstorti Exp $
source("data.m.tmp");
N==M || error("not implementeed yet!!");
U=aload("smoke.state.0.tmp");
ref = 4;

u = reshape(U(:,1),N+1,N+1);
NN = ref*N;

uu(1:ref:NN+1,1:ref:NN+1) = u;

for k=2:NN
  kk = floor((k-1)/ref)+1;
  alpha = (k-(kk-1)*ref)/ref;
  uu(k,:) = (1-alpha)*uu(kk,:) + alpha * uu(kk+1,:);
endfor

for k=2:NN
  kk = floor((k-1)/ref)+1;
  alpha = (k-kk)/ref;
  uu(:,k) = (1-alpha)*uu(:,kk) + alpha * uu(:,kk+1);
endfor
