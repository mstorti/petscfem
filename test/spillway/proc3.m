source("data.m.tmp");

fs = (1:Nx+1)'*(Ny+1);
vfs=[];
k=0;
while 1
  file = ["wave.state_" int2str(k) ".tmp"];
  [info, err, msg] = stat (file);
  if err; break; endif
  U = aload(file);
  vfs = [vfs U(fs,2)];
  k = k+1;
endwhile

xfs = aload("wave.fsh.tmp");
nfs = 100;
rem(rows(xfs),nfs)==0 || error("bad number of rows in .fsh file");
yfs = reshape(xfs(:,2),nfs,rows(xfs)/nfs);
