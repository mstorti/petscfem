###key proc.m

xfs=aload("spillway.xfs.tmp");
xfsh = [xfsh xfs(:,1)];
yfsh = [yfsh xfs(:,2)];
plot(xfsh,yfsh);
pause

system(["perl -ne 'print \"$1\n\" if m/control points: " \
	"(.+)/;' nohup.log >spillway.fs_conv.tmp"]);
h = aload("spillway.fs_conv.tmp");
plot(h);
