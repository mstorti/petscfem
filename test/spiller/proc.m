###key proc.m

xfs=aload("spiller.xfs.tmp");
xfsh = [xfsh xfs(:,1)];
yfsh = [yfsh xfs(:,2)];
plot(xfsh,yfsh);
pause

system(["perl -ne 'print \"$1\n\" if m/control points: " \
	"(.+)/;' nohup.log >spiller.fs_conv.tmp"]);
h = aload("spiller.fs_conv.tmp");
plot(h);
