n=10;
Uin=0;
Uex=1;

w=zhomo([1 1.5 0 pi/4],n+1,n+1);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);
z=unsplitc(xnod);
zz=exp(z);
xnod = splitc(zz);

iin = (1:n+1)';
iex = n*(n+1)+iin;

xin=xnod(iin,1);
yin=xnod(iin,2);
rin = sqrt(xin.^2+yin.^2);
ein = [-yin./rin xin./rin];

asave("sector.nod",xnod);
asave("sector.con",icone);

fid = fopen("sector.fix","w");
for k=1:n+1
  fprintf(fid,"%d  1  %f\n",iin(k),Uin*ein(k,1));
  fprintf(fid,"%d  2  %f\n",iin(k),Uin*ein(k,2));
endfor

for k=1:n+1
  fprintf(fid,"%d  1  %f\n",iex(k),Uex*ein(k,1));
  fprintf(fid,"%d  2  %f\n",iex(k),Uex*ein(k,2));
endfor
fclose(fid);

nnod = rows(xnod);
inlet = (1:n+1:nnod)';
outlet = inlet+n;

fid = fopen("sector.peri","w");
theta=pi/4;
coss= cos(theta);
sinn= sin(theta);
for k=1:n+1
  fprintf(fid,"-1.  %d  1     %f  %d 1  %f %d 2\n", \
          outlet(k),coss,inlet(k),-sinn,inlet(k));
  fprintf(fid,"-1.  %d  2     %f  %d 1  %f %d 2\n", \
          outlet(k),sinn,inlet(k),coss,inlet(k));
endfor

fclose(fid);
