##__INSERT_LICENSE__
## $Id: mkcouette.m,v 1.4 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

w=zhomo([0 1.5 0 1],Nx+1,Ny+1);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

usar_y=0;
if usar_y
  xnod=xnod(:,[2 1]);
  icone = icone(:,[1 4 3 2]);
endif

asave("couette.nod.tmp",xnod);
asave("couette.con.tmp",icone);

U=1;
if !usar_y
#  f=[1.1 1. 1.2  3.4 0. 0. 0.1 0.1];
  f=[1.1 0.9 1.2  3.4 0.1 0.2 0.1 0.1];
else
  f=[0. 1. 0. U  0. 0. 0.1 0.1];
endif

nnod=rows(xnod);

fixa = [];
#for k=[2 5 6 7 8]
## for k=[7 8]
for k=[]
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
## uini(:,3) = x.*(1-x)*0.1;
asave("couette.ini.tmp",uini);

fid = fopen("couette.peri.tmp","w");
for k=1:2:nnod
  for dof = 1:8
    fprintf(fid,"%f %d %d %f %d %d\n",-1.,k+1,dof,1.,k,dof);
  end
end
fclose(fid);
