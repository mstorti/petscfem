source("data.m.tmp");

if 0
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
endif

xfs = aload("wave.fsh.tmp");
nfs = Nx;
rem(rows(xfs),nfs)==0 || error("bad number of rows in .fsh file");
yfs = reshape(xfs(:,2),nfs,rows(xfs)/nfs);
xfs = reshape(xfs(:,1),nfs,rows(xfs)/nfs);

axis([min(min(xfs)) max(max(xfs)) min(min(yfs)) max(max(yfs))])
for k=1:columns(yfs)
  plot(xfs(:,k),yfs(:,k));
  pause(0.2);
endfor
