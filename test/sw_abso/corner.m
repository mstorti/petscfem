clear;
source("corner.data.m.tmp");

x=(-L:ds:0.0);
p=length(x);
xnod = [x' x' 0.*ones(p,1)];
x=(ds:ds:L);
p=length(x);
xnod = [xnod; x' 0.*ones(p,1) 0.*ones(p,1)];
nnod = size(xnod,1);

fid=fopen("corner.nod.tmp","w");
for j=1:nnod
  fprintf(fid,"%f %f %f\n",xnod(j,1),xnod(j,2),xnod(j,3));
endfor
fclose(fid);
icone1=[2:nnod-2]';
icone2=[3:nnod-1]';
icone=[icone1 icone2];
fid=fopen("corner.con.tmp","w");
for k=1:nnod-3
  fprintf(fid,"%d %d\n",icone(k,1),icone(k,2));
endfor
fclose(fid);

h_ini=h0*ones(nnod,1);
h_ini(1)=0.;
h_ini(nnod)=0.;
u_ini=u0*ones(nnod,1);
u_ini(1)=0.;
u_ini(nnod)=0.;
fid=fopen("corner.ini.tmp","w");
for i=1:nnod
fprintf(fid,"%f %f\n",u_ini(i),h_ini(i));
endfor
fclose(fid);

#rest=[1 2 3 4 nnodesn-1 u0 h0 -1.];#salida
rest = [2 3 4 5 1 u0+0.1*u0 h0+0.1*h0 -1.];
rest=[rest; nnod-1 nnod-2 nnod-3 nnod-4 nnod u0 h0 1.];#salida
fid=fopen("corner.rest.tmp","w");
fprintf(fid,"%d %d %d %d %d %f %f %f\n",rest');
fclose(fid);

if (0)
fixa=[];
fixa = [fixa; 1 1 u0;nnod-1 2 h0 ];
fid=fopen("spiral.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);
endif
clear all;