L=5;
N=50;

h=L/N;
w=zhomo([-L/2 L/2  0 h],N+1,2,[3 1 3 1 0 1]);
[xnod,icone]=pfcm2fem(w);

icone = icone(:,[1 4 3 2]);
nnod = rows(xnod);

x=xnod(:,1);
sigma=.8;
Dh=.3;
U = [.3 0 1];
U = U(ones(2*(N+1),1),:);

dh = Dh*exp(-(xnod(:,1)/sigma).^2);

xnod = [xnod dh];

asave("channel_stretch.nod",xnod);
