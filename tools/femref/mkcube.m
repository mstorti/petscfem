N=10;

w=zhomo([-1 1 -1 1],N+1,N+1);
[x2,i2] = pfcm2fem(w);

[xnod,icone] = extrude(x2,i2,N,1/N);
xnod(:,3) = 2*xnod(:,3)-1;

asave("cube.con-hexa.tmp",icone);
asave("cube.nod.tmp",xnod);

system("../hexasplit.bin -i cube.con-hexa.tmp -o cube.con.tmp");
iconet = aload("cube.con.tmp");
asave("cube.con.tmp",iconet-1);

nnod = rows(xnod);
v = pvec(xnod,[0,0,1]);
rho = l2(v);
v = leftscal(1./(0.05+rho).^2,v);

asave("cube.state.tmp",[v,zeros(nnod,1)]);

