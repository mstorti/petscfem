###key fixstate.m
### $Id: fixstate.m,v 1.3 2005/02/07 16:04:34 mstorti Exp $

file = getenv("state_file");
## file = "comp_corner_Ma_10_Euler.state.tmp";
u = aload(file);
indx = find(u(:,1)==0);
uref = u(1,:);
u(indx,:) = uref(ones(length(indx),1),:);

indx = find(isnan(u(:,1)));
ubad = [0.001,0,0,100];
u(indx,:) = ubad(ones(length(indx),1),:);

asave(file,u);
