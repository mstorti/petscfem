source("data.m.tmp");
axis_nodes = 1+(0:Nz)'*(Nr+1);
xnod = aload("vtube.nod.tmp");
z = xnod(axis_nodes,3);
if 0
  k=0;
  rho = [];
  p = [];
  w = [];
  v = [];
else
  k=columns(rho);
endif

while 1
  file = ["vtube.state_" int2str(k) ".tmp"];
  [info, erro, msg] = stat (file);
  if erro; break; endif
  printf("loading %s\n",file);
  U=aload(file);
  rho = [rho U(axis_nodes,1)];
  v = [v U(axis_nodes,3)];
  w = [w U(axis_nodes,4)];
  p = [p U(axis_nodes,5)];
  k = k+1;
endwhile
