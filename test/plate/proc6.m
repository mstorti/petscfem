##__INSERT_LICENSE__
## $Id: proc6.m,v 1.3 2003/01/09 13:39:57 mstorti Exp $
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
T=mean(TT);
St = 2/T;
maxerrel = max(abs(T-TT))/T;
printf("Strouhal %f, max deviation (rel.) %f\n",St,maxerrel);
