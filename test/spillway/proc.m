###key proc.m

#  xfs=aload("spillway.xfs.tmp");
#  xfsh = [xfsh xfs(:,1)];
#  yfsh = [yfsh xfs(:,2)];
#  plot(xfsh,yfsh);
#  pause

system(["perl -ne 'print \"$1\n\" if m/control points: " \
	"(.+)/;' nohup.log >spillway.fs_conv.tmp"]);
h = aload("spillway.fs_conv.tmp");
plot(h);

k=1;
xfsh = [];
yfsh = [];
while 1
  file = ["spillway.xfs_" int2str(k) ".tmp"];
  [info, err, msg] = stat (file);
  if err; break; endif
  xfs = aload(file);
  xfsh = [xfsh xfs(:,1)];
  yfsh = [yfsh xfs(:,2)];
  k = k+1;
endwhile
