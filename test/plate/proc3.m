##__INSERT_LICENSE__
## $Id: proc3.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
xnod = aload("cylin.nod.tmp");
nodes_x = find(abs(xnod(:,2))<1e-6 & xnod(:,1)>0);
[x,indx] = sort(xnod(nodes_x,1));
nodes_x = nodes_x(indx);
clear indx

u = [];
v = [];
p = [];
for k=0:81
  eval(["U = aload(\"qq/outvector_" num2str(k) "\");"]);
  u = [u U(nodes_x,1)];
  v = [v U(nodes_x,2)];
  p = [p U(nodes_x,3)];
endfor
clear U
