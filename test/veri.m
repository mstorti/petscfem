##__INSERT_LICENSE__
## $Id: veri.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
nn=[64 128 256 512 1024];

e_05=[];
e_1=[];
NN=nn(length(nn));
for k=nn(1:length(nn)-1)
  eval(['e_05=[e_05;max(abs(u05_' int2str(k) '-u05_' int2str(NN) '))];']); 
  eval(['e_1=[e_1;max(abs(u1_' int2str(k) '-u1_' int2str(NN) '))];']); 
endfor

e_1=e_1(:,2);
e_05=e_05(:,2);

loglog(1./nn(1:length(nn)-1),[e_1 e_05])
