u=aload("gfshock3d.dx-state-3d.tmp");
nnod = rows(u);
rem(nnod,2)==0 || error("not even number of nodes");

Rgas=287*29/27.5;
Tin = 4170;
pin = 6e5;
gamma=1.17;

cin=sqrt(1.17*Rgas*Tin);
rhoin = pin/(Rgas*Tin);

uu = rightscal(u(1:nnod/2,[1 4 2 5]),[rhoin cin cin rhoin*cin^2]);

asave("gfshock2d.dx-state.tmp",uu);
