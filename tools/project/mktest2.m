N = 10;
L = 1;

w = zhomo([0 L 0 L],N+1,N+1);
[xnod,icone] = pfcm2fem(w);

icone = [icone(:,[2 4 1]);
	 icone(:,[4 2 3])];

nnod = rows(xnod);
x = 2*xnod(:,1)-L;
z = sqrt(L^2-x.^2);
xnod = [xnod,z];

asave("mesh3.nod",xnod);
asave("mesh3.con",icone);
