N = 50;
w = zhomo([0 1 0 1],N+1,N+1);

[xnod,icone] = pfcm2fem(w);

icone = [icone(:,[2 4 1]);
	 icone(:,[4 2 3])];

asave("square1.nod",xnod);
asave("square1.con",icone);

x = xnod(:,1);
y = xnod(:,2);
phi = x.*(1-x).*y.*(1-y);

u = [phi,sqrt(phi)];
asave("square1.dat",u);

xnod2 = rand(40,2);
asave("square2.nod",xnod2);

system("./project2.bin");

xx = xnod2(:,1);
yy = xnod2(:,2);

phi = xx.*(1-xx).*yy.*(1-yy);
uana = [phi,sqrt(phi)];

uinterp = aload("pinterp.tmp");

printf("||uana|| %f, ||uinterp|| %f, ||uana-uinterp|| %f\n", merr(uana),merr(uinterp),merr(uana-uinterp));
