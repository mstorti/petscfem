n=16;

h=1/n;
w=zhomo([0 1 0 1],n+1,n+1);
[xnod,icone]=pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

N=(n+1)^2;
peri=[(1:n+1)' (n+1)*n+(1:n+1)';  # periodic in x
      (1:n+1:N)' (1:n+1:N)'+n];
nperi=rows(peri);

O=ones(nperi,1) ;
constr = [-O peri(:,1) O          +O peri(:,2) O; 
          -O peri(:,1) 2*O        +O peri(:,2) 2*O; 
          -O peri(:,1) 3*O        +O peri(:,2) 3*O];

u = [0. 0. 0.];
u=u(ones(N,1),:);

asave("homof.nod.tmp",xnod);
asave("homof.con.tmp",icone);
asave("homof.peri.tmp",constr);
asave("homof.ini.tmp",u);


