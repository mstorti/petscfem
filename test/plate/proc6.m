##__INSERT_LICENSE__
## $Id: proc6.m,v 1.5 2003/01/11 15:24:08 mstorti Exp $
TT=[];
nx = size(v,1);
for kk=min(nx,10):nx
  vv=v(kk,:)';
  nn = size(v,2);
  indx = find(vv(1:nn-1).*vv(2:nn)<0);
  if size(indx)==0; continue; endif
  indx = indx-vv(indx)./(vv(indx+1)-vv(indx));
  TT = [TT;
       2*Dt*[xcent(indx) diff(indx)]];
endfor
if rows(TT)==0 
  disp("Couldn't find any periods to make statistics");
  return
endif
[bid,indx] = sort(TT(:,1));
clear bid
TT=TT(indx,:);
t=TT(:,1);
TT=TT(:,2);
T=mean(TT);
St = 2/T;
maxerrel = max(abs(T-TT))/T;
printf("Strouhal %f, max deviation (rel.) %f\n",St,maxerrel);
