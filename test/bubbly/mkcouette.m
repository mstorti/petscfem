source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

w=zhomo([0 1 0 1],N+1,2);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

asave("couette.nod.tmp",xnod);
asave("couette.con.tmp",icone);

f=[0. 1. 0. 0. 0. 0. 0.1 0.1];

nnod=rows(xnod);

fixa = [];
for k=[2 5 6 7 8]
  fixa = [fixa;
	  (1:nnod)' k*ones(nnod,1) f(k)*ones(nnod,1)];
end

fid = fopen("couette.fixa.tmp","w");
for kk=1:rows(fixa);
  fprintf(fid,"%d %d %f\n",fixa(kk,1),fixa(kk,2),fixa(kk,3));
end
fclose(fid);

uini = kron(ones(nnod,1),f);
x = xnod(:,1);
uini(:,3) = x.*(1-x)*0.1;
asave("couette.ini.tmp",uini);

fid = fopen("couette.peri.tmp","w");
for k=1:2:nnod
  for dof = 1:8
    fprintf(fid,"%f %d %d %f %d %d\n",-1.,k+1,dof,1.,k,dof);
  end
end
fclose(fid);
