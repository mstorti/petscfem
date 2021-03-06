source("data.m.tmp");
h = 2/N;
dd = 0.1;

w=zhomo([-1 1 -1 1],N+1,N+1);
[x2,i2] = pfcm2fem(w);
if 0
  alpha = pi/4;
  R = [cos(alpha), -sin(alpha);
       sin(alpha), cos(alpha)];
  x2 = x2*R';
endif

[xnod,icone] = extrude(x2,i2,N,1/N);
xnod(:,3) = 2*xnod(:,3)-1;

asave("cube.con-hexa.tmp",icone);
asave("cube.nod.tmp",xnod);

system("../hexasplit.bin -i cube.con-hexa.tmp -o cube.con.tmp");

icone = aload("cube.con.tmp");
asave("cube.con0.tmp",icone-1);

nnod = rows(xnod);
if 1
  v = pvec(xnod,[0,0,1]);
  rho = l2(v);
  v = leftscal(1./(dd^2+rho.^2),v);
else
  v = pvec([0,0,1],xnod);
endif
asave("cube.state.tmp",[v,zeros(nnod,1)]);

nelem = rows(icone);
asave("cube.state-elem.tmp",5*ones(nelem,3));

return

system("make getsurf");

grad_ue = aload("cube.grad-u.tmp");

we = [+grad_ue(:,10)-grad_ue(:,7), \
     -grad_ue(:,3)+grad_ue(:,9), \
     +grad_ue(:,5)-grad_ue(:,2)];
awe = l2(we);

asave("cube.awe.tmp",awe);

grad_u = aload("cube.grad-un.tmp");

w = [+grad_u(:,10)-grad_u(:,7), \
     -grad_u(:,3)+grad_u(:,9), \
     +grad_u(:,5)-grad_u(:,2)];
aw = l2(w);

asave("cube.aw.tmp",aw);

if 0
  xe = pfnd2ele(xnod,surf_con,xnod);

  rhoe=l2(xe(:,1:2));
  awt = 2*dd^2./(dd^2+rhoe.^2).^2;
  dudr = (dd^2-rhoe.^2)./(dd^2+rhoe.^2).^2;

  ## indx = find(xe(:,2)<h & xe(:,2)>0 & xe(:,3)>0 & abs(xe(:,1))<0.99);
  indx = find(abs(xe(:,2))<1e-4 & xe(:,3)>0);
  indx=sortby(xe(indx,1),indx);
  plot(xe(indx,1),[w(indx,3),awt(indx)])
endif

