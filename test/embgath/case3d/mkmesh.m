## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");

h = L/N;
w=zhomo([0 h 0 L],2,N+1,[1 0 1 1 hratio 1]);

[x2,i2]=pfcm2fem(w);
i2 = i2(:,[1,4,3,2]);

[x3d,ico3d] = extrude(x2,i2,1,h);

asave("strip3d.nod.tmp",x3d);
asave("strip3d.con.tmp",ico3d);

tol=1e-5;
bot = (N+1)*(0:3)'+1;
top = bot+N;

pffixa("strip3d.fixa-wall.tmp",[bot;top],1:3);
i1 = (1:N+1)';
i2 = i1+(N+1);
i3 = i2+(N+1);
i4 = i3+(N+1);
pfperi("strip3d.peri.tmp",i2,i1,(1:4));
pfperi("strip3d.peri.tmp",i3,i1,(1:4),"a");
pfperi("strip3d.peri.tmp",i4,i1,(1:4),"a");

## Connectivity
icosurf = [([1,2,4,3]-1)*(N+1)+1;
           [1,3,4,2]*(N+1)];
if !use_exterior_normal;
  icosurf = icosurf(:,[1,4,3,2]);
endif
asave("strip3d.surf-con.tmp",[icosurf,ones(2,8)]);
