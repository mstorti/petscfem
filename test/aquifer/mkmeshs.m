source("data.m.tmp");

x=[(0:Nx)'/Nx*Lx zeros(Nx+1,1)];
ico=[(1:Nx)' (2:Nx+1)'];

asave("stream.nod.tmp",x);
asave("stream.con.tmp",ico);
