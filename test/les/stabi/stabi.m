##__INSERT_LICENSE__
## $Id: stabi.m,v 1.4 2003/01/08 15:49:04 mstorti Exp $
#source("ini.data");
#system("make run");

petsc_data_name="system.dat";
petscload;

h=Ly/N;
n=(N+1)*3-4;
M10=getblock(atet,n+(1:n)',(1:n)',"f");
M11=getblock(atet,n+(1:n)',n+(1:n)',"f");
M12=getblock(atet,n+(1:n)',2*n+(1:n)',"f");

petsc_data_name="flux_jaco.dat";
petscload;

C10=getblock(atet,n+(1:n)',(1:n)',"f");
C11=getblock(atet,n+(1:n)',n+(1:n)',"f");
C12=getblock(atet,n+(1:n)',2*n+(1:n)',"f");

M10=M10-C10;
M11=M11-C11;
M12=M12-C12;

ip=[1:3:n-1 n]';
iu=[2:3:n-1]';
iv=[3:3:n-1]';
iuv=[iu;iv];
Nu=length(iu);
Nv=length(iv);

lambda1=1;
lambda2=10;
nla=50;
if 1
  lav=logspace(log10(lambda1),log10(lambda2),nla)';
else
  nla=1;
  lav=200;
endif
alpha=zeros(nla,1);
for k=1:nla
  kw=2*pi/lav(k);
  phase = exp(i*kw*h);
  M=phase*M10+M11+M12/phase;
  C=phase*C10+C11+C12/phase;
  
  Muu=M(iuv,iuv);
  Mpu=M(ip,iuv);
  K=-C(iuv,iuv);
  G=-C(iuv,ip);
  H=-C(ip,iuv);
  D=-C(ip,ip);

  MM=Muu-G*(D\Mpu);
  KK=K-G*(D\H);
  
  d=eig(MM\KK);
  [bid,j]=max(real(d));
  alpha(k) = d(j);
endfor

if !stokes  # general case
  [v,d]=eig(MM\KK);
  d=diag(d);
  [bd,indx]=sort(-real(d));
  d=d(indx);
  v=v(:,indx);
  ph=(v(1,:)./abs(v(1,:))).';
  v=rightscal(v,1./ph);
  clear ph
  uu=v(1:Nu,:);
  vv=-i*v(Nu+(1:Nv),:);
else
  [v,d]=eig(MM\KK);
  d=takereal(diag(d));
  [d,indx]=sort(-d);
  v=v(:,indx);
  ph=(v(1,:)./abs(v(1,:))).';
  v=rightscal(v,1./ph);
  clear ph
  uu=takereal(v(1:Nu,:));
  vv=takereal(-i*v(Nu+(1:Nv),:));
endif

return

[NM,R]=nullrange(M');
NM=NM';

RM=R*M;
RC=R*C;
NC=NM*C;

return

## this can serve to check if the periodic boundary conditions are
## correctly set (should be C00=C11=C22, C10=C21=C02, etc...)

#  C00=getblock(atet,(1:n)',(1:n)',"f");
#  C01=getblock(atet,(1:n)',n+(1:n)',"f");
#  C02=getblock(atet,(1:n)',2*n+(1:n)',"f");

#  C22=getblock(atet,2*n+(1:n)',2*n+(1:n)',"f");

#for signo=[-1 1]
for signo=1
  phase = exp(signo*i*k*h);
  [v,d]=eig(phase*C10+C11+(1/phase)*C12);
  d=diag(d);
  [d,indx]=sort(real(d));
  v=v(:,indx);
  min(real(d))
  plot([real(d) imag(d)])
  if signo==-1
    pause
  endif
endfor
return

indx=find(abs(d)>1e-3);
d=d(indx);
v=v(:,indx);

[d,indx]=sort(d);
v=v(:,indx);
