##__INSERT_LICENSE__
## $Id: mmesh.m,v 1.5 2003/11/16 13:18:50 mstorti Exp $
source("data.m.tmp");

w=zhomo([0 1 0 1/N],N+1,2);
[xnod,icone]=pfcm2fem(w);

icone=icone(:,[1 4 3 2]);
asave("burgers.nod.tmp",xnod);
asave("burgers.con.tmp",icone);

x=xnod(:,1);
u=1-2*x;
asave("burgers.ini.tmp",u);

u=ones(size(u));
asave("burgers.inip.tmp",u);
asave("burgers.inim.tmp",u);
