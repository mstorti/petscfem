###key fixstate.m
### $Id: fixstate.m,v 1.1 2005/02/02 22:56:04 mstorti Exp $

u = aload("cylabso.dx-state.tmp");
indx = find(u(:,1)==0);
u(indx,[1,4]) = 1e-5;
asave("cylabso.dx-state.tmp",u);
