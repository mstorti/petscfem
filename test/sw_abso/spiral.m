#reference nodes:1 and nnod
clear;
source("spiral.data.m.tmp");
t = (0:nt:k*2*pi);
xnod = [R*cos(t)' R*sin(t)' -1*t'];
pp=xnod(rows(xnod),3);
nnod = size(xnod,1);

fid=fopen("spiral.nod.tmp","w");
for j=1:nnod
  fprintf(fid,"%f %f %f\n",xnod(j,1),xnod(j,2),xnod(j,3));
endfor
#fprintf(fid,"%f %f %f\n",xnod(1,1), -nt, xnod(1,3));
#fprintf(fid,"%f %f %f\n", xnod(1,1), nt , pp);
fclose(fid);
icone1=[2:nnod-2]';
icone2=[3:nnod-1]';
icone=[icone1 icone2];
fid=fopen("spiral.con.tmp","w");
for k=1:nnod-3
  fprintf(fid,"%d %d\n",icone(k,1),icone(k,2));
endfor
fclose(fid);

h_ini=h0*ones(nnod,1);
h_ini(2)=h0*1.1;
h_ini(3)=h0*1.07;
h_ini(4)=h0*1.05;
h_ini(5)=h0*1.02;
h_ini(1)=0.;
h_ini(nnod)=0.;
u_ini=u0*ones(nnod,1);
u_ini(1)=0.;
u_ini(nnod)=0.;
fid=fopen("spiral.ini.tmp","w");
for i=1:nnod
fprintf(fid,"%f %f\n",u_ini(i),h_ini(i));
endfor
fclose(fid);

rest = [3 4 5 6 1 2 -1.];
rest=[rest; nnod-2 nnod-3 nnod-4 nnod-5 nnod nnod-1 1.];#salida
fid=fopen("spiral.rest.tmp","w");
fprintf(fid,"%d %d %d %d %d %d %f\n",rest');
fclose(fid);

if (1)
fixa=[];
fixa = [fixa; 2 1 u0;2 2 h0; nnod-1 1 u0; nnod-1 2 h0 ];
fid=fopen("spiral.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);
endif
clear all;