##$Id: proc3.m,v 1.4 2003/06/08 13:10:46 mstorti Exp $
source("data.m.tmp");
N = M = 160;
ref = 2;
NN = ref*N;

if 0
  V = aload("sqcav.state.N80");
  VV = zeros((NN+1)^2,1);
  for k=1:3
    VV(:,k) = refin(V(:,k),ref);
  endfor
  asave("sqcav.state.N160",VV);
endif

N==M || error("not implemented yet!!");
for s=0:31
  s
  U=aload(["STEPS/smoke.state." int2str(s) ".tmp"]);
  UU = zeros((NN+1)^2,1);
  UU = refin(U,ref);
  asave(["STEPS_REF/smoke.state." int2str(s) ".tmp"],UU);
endfor
