system("./cpustat.pl nohup.log > stat.dat");
a=aload("stat.dat");
a = [a sum(a,2)];
a(1,:) = [];
t = mean(a);
n = rows(a);
dt = sqrt(sum(a.^2)/n-t.^2);
[t' dt']
