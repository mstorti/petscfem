nnod_axis = 101;
V = aload("cylin.axis_nodes.rslt.tmp");

rem(rows(V),nnod_axis)==0 || \
    error("not correct number of rows in cylin.axis_nodes.rslt.tmp");

nsteps = round(rows(V)/nnod_axis);
u = reshape(V(:,2),nnod_axis,nsteps);
v = reshape(V(:,3),nnod_axis,nsteps);
p = reshape(V(:,4),nnod_axis,nsteps);
clear V
