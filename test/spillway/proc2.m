fs = aload("spillway.nod_fs.tmp");
xnod = aload("spillway.nod.tmp");
xfs = xnod(fs,1);
pfs = [];
for k=0:5
  state = aload(sprintf("spillway.state_%d.tmp",k));
  pfs = [pfs state(fs,3)];
endfor
