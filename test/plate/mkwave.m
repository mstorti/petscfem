##__INSERT_LICENSE__
## $Id: mkwave.m,v 1.2 2003/02/24 23:19:34 mstorti Exp $
global H L 

source("data.m.tmp");
Hy = H;

XNOD = [1 0 0;
	2 L 0;
	3 L H;
	4 0 H];

XNOD = XNOD(:,2:3);

ICONE = [1 2 3 4];

H = [1 2 Nx 1 0 xratio;
     2 3 Ny 1 0 1];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H);
inlet = mesher_bound(mesh,[1 4]);
outlet = mesher_bound(mesh,[2 3]);
bottom = mesher_bound(mesh,[1 2]);
top = mesher_bound(mesh,[4 3]);

asave("wave.nod.tmp",xnod);
asave("wave.con.tmp",icone);

fid = fopen("wave.fixa_in.tmp","w");
for k=inlet'
  uu = 1+du*cos(pi*xnod(k,2)/Hy);
  fprintf(fid,"%d %d  %f\n",k,1,uu);
endfor
fclose(fid);

x = xnod((0:2)'*(Ny+1)+1,1);
cx = cloud(x,1,2);
hy = Hy/Ny;

fid2 = fopen("wave.null_w.tmp","w");
for k=2:Ny
  k1 = k  + (Ny+1);
  k2 = k1 + (Ny+1);
  fprintf(fid2,"   %f   %d %d",1/(2*hy),k+1,1);
  fprintf(fid2,"   %f   %d %d",-1/(2*hy),k-1,1);
  fprintf(fid2,"   %f   %d %d",-cx(1),k,2);
  fprintf(fid2,"   %f   %d %d",-cx(2),k1,2);
  fprintf(fid2,"   %f   %d %d",-cx(3),k2,2);
  fprintf(fid2,"\n");
endfor
fclose(fid2);

fid = fopen("wave.fixa_out.tmp","w");
for k=outlet'
  fprintf(fid,"%d %d  %f\n",k,2,0.);
endfor
fclose(fid);

fid = fopen("wave.fixa_top.tmp","w");
for k=top'
  fprintf(fid,"%d %d  %f\n",k,2,0.);
endfor
fclose(fid);

fid = fopen("wave.fixa_bot.tmp","w");
for k=bottom'
  fprintf(fid,"%d %d  %f\n",k,2,0.);
endfor
fclose(fid);

uini = [1 0 0];
uini = uini(ones(rows(xnod),1),:);
asave("wave.ini.tmp",uini);

fid = fopen("wave.coupler.tmp","w");
for k=(length(inlet)-1):-1:1
  fprintf(fid,"%d %d\n",inlet(k+1),inlet(k));
endfor
fclose(fid);
