##$Id: proc3.m,v 1.2 2003/06/06 18:12:52 mstorti Exp $
source("data.m.tmp");
ref = 4;

N==M || error("not implementeed yet!!");
for s=0:0

  U=aload(["STEPS/smoke.state." int2str(s) ".tmp"]);
  NN = ref*N;
  UU = zeros((NN+1)^2);
  for k=1:4
    UU(:,k) = refin(U(:,k),ref);
  endfor
  asave(["STEPS_REF/smoke.state." int2str(s) ".tmp"]);

endfor

