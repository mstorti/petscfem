output
d=max(abs(D));
Iu = find(d==0);
Ip = find(d>0);

clear d
disp('A:'); A=takef(A,Iu,Iu);
disp('B:'); B=takef(B,Iu,Ip);
disp('C:'); C=takef(C,Ip,Iu);
disp('D:'); D=takef(D,Ip,Ip);

disp('resp:'); resp=takef(resp,Ip);
disp('resu:'); resu=takef(resu,Iu);

nu=length(Iu);
np=length(Ip);

u=zeros(nu,1);
p=zeros(np,1);
omega=0.5;
h=[];
while 1
  unew=A\(resu-B*p);
  pnew=D\(resp-C*u);

  du=norm(unew-u);
  dp=norm(pnew-p);
  disp(sprintf('du = %f, dp = %f, du+dp = %f',du,dp,du+dp));
  h=[h;
     du dp du+dp];

  u=u+omega*(unew-u);
  p=p+omega*(pnew-p);
end
