##__INSERT_LICENSE__
## $Id: proc5.m,v 1.12 2003/03/01 23:26:19 mstorti Exp $
source("data.m.tmp");
tmp = aload("cylin.axis_nodes.tmp");
xnod = aload("cylin.nod.tmp");
x = xnod(tmp,1);
nnod_axis = length(tmp);
clear tmp
V = aload("cylin.axis_nodes.rslt.tmp");

rem(rows(V),nnod_axis)==0 || \
    error("not correct number of rows in cylin.axis_nodes.rslt.tmp");

nsteps = round(rows(V)/nnod_axis);
u = reshape(V(:,2),nnod_axis,nsteps);
v = reshape(V(:,3),nnod_axis,nsteps);
p = reshape(V(:,4),nnod_axis,nsteps);
clear V

if 1
#  v(:,1:240)=[];
  v = v(:,500:1100);
  v = v(50:110,:);
  proc6
endif
