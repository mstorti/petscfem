n=20;

h=1/n;
w=zhomo([0 1 0 1],n+1,n+1,[1 5 1 1 5 1]);
[xnod,icone]=pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

N=(n+1)^2;

fixa0 = [(1:n+1), \
         (1:n+1:N), \
         n*(n+1)+1:1:N]';

fixa1 = [2*(n+1):n+1:n*(n+1)]';

fixag(fixa0,1,0,"sqcav.fixau0.tmp");
fixag(fixa1,1,1,"sqcav.fixau1.tmp");
fixag([fixa0;fixa1],2,0,"sqcav.fixav0.tmp");

#  O=ones(nperi,1) ;
#  constr = [-O peri(:,1) O          +O peri(:,2) O; 
#            -O peri(:,1) 2*O        +O peri(:,2) 2*O; 
#            -O peri(:,1) 3*O        +O peri(:,2) 3*O];

#  u = [0. 0. 0.];
#  u=u(ones(N,1),:);

asave("sqcav.nod.tmp",xnod);
asave("sqcav.con.tmp",icone);
