##__INSERT_LICENSE__
## $Id: proc4.m,v 1.4 2003/01/08 15:49:04 mstorti Exp $
if 1
  x = aload("cylin.nod.tmp");
  nod_ext = aload("ext.coupling_nodes.tmp");
  nod_ext = nod_ext(:,1);
  x = x(nod_ext,1);
endif
  
if !exist("U");
  U = [];
endif

step=columns(U);
while 1
  file = ["cylin.state_" int2str(step) ".tmp"];
  [info, err, msg] = stat (file);
  if err; break; endif
  printf("loading step %d\n",step);
  UU = aload(file);
  U = [U UU(nod_ext,1)];
  step = step + 1;
endwhile
clear UU file
