step=columns(U);
nod_ext = aload("ext.coupling_nodes.tmp");
nod_ext = nod_ext(:,1);
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
