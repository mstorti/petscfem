##__INSERT_LICENSE__
## $Id: proc6.m,v 1.2.2.1 2003/01/11 13:56:26 mstorti Exp $
TT=[];
nx = size(v,1);
for kk=min(nx,10):nx
  vv=v(kk,:)';
  nn = size(v,2);
  indx = find(vv(1:nn-1).*vv(2:nn)<0);
  if size(indx)==0; continue; endif
  indx = indx-vv(indx)./(vv(indx+1)-vv(indx));
  TT = [TT;
       2*Dt*diff(indx)];
endfor
if length(TT)==0 
  disp("Couldn't find any periods to make statistics");
  return
endif

T=mean(TT);
St = 2/T;
maxerrel = max(abs(T-TT))/T;
printf("Strouhal %f, max dev. (rel.) %f\n",St,maxerrel);
