##__INSERT_LICENSE__
## $Id: proc5.m,v 1.9 2003/01/16 15:31:54 mstorti Exp $
source("data.m.tmp");
tmp = aload("cylin.axis_nodes.tmp");
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
  v(:,1:300)=[];
  v = v(50:70,:);
  proc6
endif
