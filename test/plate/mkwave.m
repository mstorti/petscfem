##__INSERT_LICENSE__
## $Id: mkwave.m,v 1.1 2003/02/24 21:06:12 mstorti Exp $
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

fid_pf = fopen("wave.press_filt.tmp","w");
fid = fopen("wave.fixa_in.tmp","w");
for k=inlet'
  uu = 1+du*cos(pi*xnod(k,2)/Hy);
  p = -0.5*uu^2;
  k1 = k  + (Ny+1);
  k2 = k1 + (Ny+1);
  fprintf(fid,"%d %d  %f\n",k,1,uu);
  fprintf(fid,"%d %d  %f\n",k,3,p);
  fprintf(fid_pf,"%f  %d %d    %f  %d %d    %f  %d %d\n",
	  -1.,k,3,2,k1,3,-1,k2,3);
endfor
fclose(fid);
fclose(fid_pf);

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
