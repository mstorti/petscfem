system("./cpustat.pl nohup-fs-N125.log > stat.dat.tmp");
a=aload("stat.dat.tmp");
unlink("stat.dat.tmp");
a = [a sum(a,2)];
a(1,:) = [];
t = mean(a);
n = rows(a);
dt = sqrt(sum(a.^2)/n-t.^2);
[t' dt']
