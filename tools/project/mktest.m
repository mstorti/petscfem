N = 100;
nnod2 = 4000;
w = zhomo([0 1 0 1],N+1,N+1);

[xnod,icone] = pfcm2fem(w);

icone = [icone(:,[2 4 1]);
	 icone(:,[4 2 3])];

asave("square1.nod.tmp",xnod);
asave("square1.con.tmp",icone);

x = xnod(:,1);
y = xnod(:,2);
phi = x.*(1-x).*y.*(1-y);

u = [phi,sqrt(phi)];
asave("square1.dat.tmp",u);

xnod2 = rand(nnod2,2);
asave("square2.nod.tmp",xnod2);
return

start = clock();
system("./project2.bin");
printf("elapsed %f\n",etime(clock(),start));

xx = xnod2(:,1);
yy = xnod2(:,2);

phi = xx.*(1-xx).*yy.*(1-yy);
uana = [phi,sqrt(phi)];

uinterp = aload("pinterp.tmp");

printf("||uana|| %f, ||uinterp|| %f, ||uana-uinterp|| %f\n", merr(uana),merr(uinterp),merr(uana-uinterp));
