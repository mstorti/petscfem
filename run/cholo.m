N=20,
L=1.5;
DH=.2;
R=1;
tol=1e-4;

w=zhomo([-L L 0 L],2*N+1,N+1);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

x = xnod(:,1);
y = xnod(:,2);
r= sqrt(x.^2+y.^2);
H = DH*(1-(r/R).^2);
H = H.*(H>0);

xnod = [xnod H];

indx = [find(abs(xnod(:,2))<tol);
	find(abs(xnod(:,2)-L)<tol);];

asave("cholo.nod",xnod);
asave("cholo.con",icone);

fid = fopen("cholo.fixa","w");
for k=1:length(indx)
  fprintf(fid,"%d %d %f \n",indx(k),2,0.);
endfor
fclose(fid);

fid = fopen("cholo.abso","w");

indx = find(abs(xnod(:,1)-L)<tol);
for k=1:length(indx)
  fprintf(fid,"%d  %f %f \n",indx(k),1.,0.);
endfor

indx = find(abs(xnod(:,1)+L)<tol);
for k=1:length(indx)
  fprintf(fid,"%d  %f %f \n",indx(k),-1.,0.);
endfor
fclose(fid);
