source("ini.data");

n=(N+1)*3-4;
C10=getblock(atet,n+(1:n)',(1:n)',"f");
C11=getblock(atet,n+(1:n)',n+(1:n)',"f");
C12=getblock(atet,n+(1:n)',2*n+(1:n)',"f");

ip=[1:3:n-1 n]';
iu=[2:3:n-1]';
iv=[3:3:n-1]';

## this can serve to check if the periodic boundary conditions are
## correctly set (should be C00=C11=C22, C10=C21=C02, etc...)

#  C00=getblock(atet,(1:n)',(1:n)',"f");
#  C01=getblock(atet,(1:n)',n+(1:n)',"f");
#  C02=getblock(atet,(1:n)',2*n+(1:n)',"f");

#  C22=getblock(atet,2*n+(1:n)',2*n+(1:n)',"f");

[v,d]=eig(C10+C11+C12);
d=diag(d);

indx=find(abs(d)>1e-3);
d=d(indx);
v=v(:,indx);

[d,indx]=sort(d);
v=v(:,indx);
