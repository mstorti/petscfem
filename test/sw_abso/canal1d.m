#$Id: canal1d.m,v 1.2.2.1 2003/10/20 15:02:59 mstorti Exp $
#crea malla y condiciones de borde e iniciales  de sw-1d
#condiciones absorbentes
## Number os elements in x, acordarse de poner nodos dummy
clear;
source("canal1d.data.m.tmp");
slope=0.0; % 10% de la long del rio
x=L*onedstr([.1 .1 .1],nx+1);# n de nodos deberia ser impar para que el bump quede en el centro
nnodes=length(x);
y=zeros(nnodes,1);
z=(-x*sqrt(2)*4e-5)*slope;
nod=[x y z];
nod=[nod;(nod(1,1)-L/(nnodes-1)) 0. 0.;(nod(nnodes,1)+L/(nnodes-1)) 0. 0.];
nnodesn=size(nod,1);
icone1=[1:nnodes-1]';
icone2=[2:nnodes]';
icone=[icone1 icone2];

fid=fopen("canal1d.nod.tmp","w");
for j=1:nnodesn
  fprintf(fid,"%f %f %f\n",nod(j,1),nod(j,2),nod(j,3));
endfor
fclose(fid);

fid=fopen("canal1d.con.tmp","w");
for k=1:nnodes-1
  fprintf(fid,"%d %d\n",icone(k,1),icone(k,2));
endfor
fclose(fid);


#parabolic bump paper mario for initial state
xc=nod(round(nx/2),1);
mu=L/2;
sig=0.04*L;
rc=.1*L;
A=3; % ahora es la altura del bump en el centro, no una relacion de aspecto.
#rmax=0.;
for i=1:nnodes,
h_ini(i)=1.+A*((1/(sqrt(2.*pi)))*exp(-0.5*((nod(i,1)-mu)/sig)^2));
endfor
h_ini(nnodesn-1)=0.;h_ini(nnodesn)=0.;
fid=fopen("canal1d.ini.tmp","w");
u0v=u0*ones(nnodes,1);
u0v(nnodesn-1)=0.;u0v(nnodesn)=0.;
for i=1:nnodesn
fprintf(fid,"%f %f\n",u0v(i),h_ini(i));
endfor
fclose(fid);

rest=[1 2 3 4 nnodesn-1 u0 h0 -1.];#salida
rest=[rest; nnodes nnodes-1 nnodes-2 nnodes-3 nnodesn u0 h0 1.];#salida
fid=fopen("canal1d.rest.tmp","w");
fprintf(fid,"%d %d %d %d %d %f %f %f\n",rest');
fclose(fid);
clear all;
