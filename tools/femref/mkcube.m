N=1;

w=zhomo([0 1 0 1],N+1,N+1);
[x2,i2] = pfcm2fem(w);

[xnod,icone] = extrude(x2,i2,N,1/N);

asave("cube.con-hexa.tmp",icone);
asave("cube.nod.tmp",xnod);

system("../hexasplit.bin -i cube.con-hexa.tmp -o cube.con.tmp");
iconet = aload("cube.con.tmp");
asave("cube.con.tmp",iconet-1);

nnod = rows(xnod);
asave("cube.state.tmp",[xnod,zeros(nnod,1)]);

