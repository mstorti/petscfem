clear;
source("centerflux.data.m.tmp");
slope=0.0; % 10% de la long del rio
x=L*onedstr([.1 .1 .1],nx+1);# n de nodos deberia ser impar para que el bump quede en el centro
nnodes=length(x);
y=zeros(nnodes,1);
z=(-x*sqrt(2)*4e-5)*slope;
nod=[x y z];
icone1=[1:nnodes-1]';
icone2=[2:nnodes]';
icone=[icone1 icone2];

fid=fopen("centerflux.nod.tmp","w");
for j=1:nnodes
  fprintf(fid,"%f %f %f\n",nod(j,1),nod(j,2),nod(j,3));
endfor
fclose(fid);

fid=fopen("centerflux.con.tmp","w");
for k=1:nnodes-1
  fprintf(fid,"%d %d\n",icone(k,1),icone(k,2));
endfor
fclose(fid);

h_inil=h0*ones(floor(nnodes/2),1);
u_inil=u0*ones(floor(nnodes/2),1);
h_inir=h0*ones(nnodes-(nnodes/2),1);
u_inir=5.*u0*ones(nnodes-(nnodes/2),1);
h_ini=[h_inil;h_inir];
u_ini=[u_inil;u_inir];
fid=fopen("centerflux.ini.tmp","w");
for i=1:nnodes
fprintf(fid,"%f %f\n",u_ini(i),h_ini(i));
endfor
fclose(fid);

fixa=[1 1 u0; 1 2 h0; nnodes 1 5*u0; nnodes 2 h0];
fid=fopen("centerflux.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);
clear all;
