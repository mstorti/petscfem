##$Id: proc3.m,v 1.5 2003/06/09 21:14:11 mstorti Exp $
source("data.m.tmp");
N = M = 80;
ref = 4;
NN = ref*N;

if 1
  V = aload("sqcav.state.N80");
  VV = zeros((NN+1)^2,1);
  for k=1:3
    VV(:,k) = refin(V(:,k),ref);
  endfor
  asave("sqcav.state.N320",VV);
endif

N==M || error("not implemented yet!!");
for s=82:113
  s
  U=aload(["STEPS/smoke.state." int2str(s) ".tmp"]);
  UU = zeros((NN+1)^2,1);
  UU = refin(U,ref);
  asave(["STEPS_REF/smoke.state." int2str(s-82) ".tmp"],UU);
endfor
