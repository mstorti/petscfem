N = 100;
L = 1;

w = zhomo([0 L 0 L],N+1,N+1);
[xnod,icone] = pfcm2fem(w);

icone = [icone(:,[2 4 1]);
	 icone(:,[4 2 3])];

nnod = rows(xnod);

if 0
  phi = 2*pi*xnod(:,1)/L;
  y = xnod(:,2);
  z = L*cos(phi);
  x = L*sin(phi);
else
  x = xnod(:,1);
  y = xnod(:,2);
  z = sqrt(L^2-x.^2);
endif
xnod = [x,y,z];

asave("mesh3.nod",xnod);
asave("mesh3.con",icone);
