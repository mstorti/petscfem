##__INSERT_LICENSE__
## $Id: mkstrip.m,v 1.1 2003/01/19 02:25:49 mstorti Exp $

source("data.m.tmp");

Dr_av = (R0-Rin)/Nr;
Dz = Dr_av;
w = zhomo([0 Dz Rin R0],2,Nr+1,[1 0 1 r_ratio 0 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);
xnod = xnod(:,[2 1]);
Dphi = Dr_av/sqrt(R0*Rin);
[x3d,ic3d] = extrude(xnod,icone,1,Dphi);
rho = x3d(:,1);
z = x3d(:,2);
phi = x3d(:,3);
x3d = [rho.*cos(phi) rho.*sin(phi) z];

asave("strip.nod.tmp",x3d);
asave("strip.con.tmp",ic3d);

fid = fopen("strip.fixa.tmp","w");
u_dof = 2;
for l=1:3, 
  fprintf(fid,"%d %d   %f\n",1,u_dof+l-1,0.); 
  fprintf(fid,"%d %d   %f\n",Nr+1,u_dof+l-1,0.); 
endfor
fclose(fid);

fid = fopen("strip.peri.tmp","w");
c = cos(Dphi); s = sin(Dphi); 
for k=1:Nr+1
  for p=1:3
    peri = k+p*(Nr+1);
    if p==1 
      for l=1:5
	fprintf(fid,"%f   %d %d   %f  %d %d\n",-1.,peri,l,1.,k,l); 
      endfor
    else
      for l=[1 4 5];
	fprintf(fid,"%f   %d %d   %f  %d %d\n",-1.,peri,l,1.,k,l); 
      endfor
      ## U_peri = cos(Dphi) *  U_k - sin(Dphi) *  V_k 
      ## V_peri = sin(Dphi) *  U_k + cos(Dphi) *  V_k 
      fprintf(fid,"%f   %d %d   %f  %d %d  %f  %d %d  \n",
	      -1.,peri,2,+c,k,2,-s,k,3); 
      fprintf(fid,"%f   %d %d   %f  %d %d  %f  %d %d  \n",
	      -1.,peri,3,+s,k,2,+c,k,3); 
    endif
  endfor
endfor
fclose(fid);

nnod = rows(x3d);
uini = ones(nnod,5)*diag([rho_in 0 0 0 p_in]);
uini(:,2:3) = [-W*x3d(:,2) +W*x3d(:,1)];
asave("strip.ini.tmp",uini);
