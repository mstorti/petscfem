N=100;

w=zhomo([0 1 0 1/N],N+1,2);
[xnod,icone]=pfcm2fem(w);

icone=icone(:,[1 4 3 2]);
asave("plano.nod",xnod);
asave("plano.con",icone);

x=xnod(:,1);
u=1-2*x;
asave("burgers.ini",u);
