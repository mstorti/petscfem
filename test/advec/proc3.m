##$Id: proc3.m,v 1.3 2003/06/06 20:43:41 mstorti Exp $
source("data.m.tmp");
N=M=80;
ref = 4;
NN = ref*N;

if 0
  V = aload("sqcav.state.N80.tmp");
  VV = zeros((NN+1)^2,1);
  for k=1:3
    VV(:,k) = refin(V(:,k),ref);
  endfor
  asave("sqcav.state.N320.tmp",VV);
endif

N==M || error("not implementeed yet!!");
for s=1:31

  s
  U=aload(["STEPS/smoke.state." int2str(s) ".tmp"]);
  UU = zeros((NN+1)^2,1);
  UU = refin(U,ref);
  asave(["STEPS_REF/smoke.state." int2str(s) ".tmp"],UU);

endfor

