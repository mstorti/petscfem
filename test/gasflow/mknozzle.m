global Rin Rn nw L

source("data.m.tmp");

XNOD = [1 -L 0;
	2 +L 0;
	3 +L Rin;
	4 -L Rin];

XNOD = XNOD(:,2:3);

ICONE = [1 2 3 4];

H = [1 2 Nx r_ratio 1 r_ratio;
     4 1 Nr 1       0 x_ratio];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H);
asave("nozzle.nod.tmp",xnod);
asave("nozzle.con.tmp",icone);

inlet = mesher_bound(mesh,[4 1]);
outlet = mesher_bound(mesh,[2 3]);
axis_b = mesher_bound(mesh,[1 2]);
axis_b = mesher_bound(mesh,[1 2]);
wall = mesher_bound(mesh,[3 4]);

fid = fopen("nozzle.fixa_inlet.tmp","w");
for k=2:length(inlet)
  node = inlet(k);
  fprintf(fid,"%d %d  %f\n",node,1,rho_ref);
  fprintf(fid,"%d %d  %f\n",node,2,u_ref);
  fprintf(fid,"%d %d  %f\n",node,3,v_ref);
endfor
fclose(fid);

fid = fopen("nozzle.fixa_wall.tmp","w");
for k=1:length(wall)
  node = wall(k);
  fprintf(fid,"%d %d  %f\n",node,2,0.);
  fprintf(fid,"%d %d  %f\n",node,3,0.);
endfor
fclose(fid);

fid = fopen("nozzle.fixa_axis.tmp","w");
for k=2:length(axis_b)-1
  node = axis_b(k);
  fprintf(fid,"%d %d  %f\n",node,3,0.);
endfor
fclose(fid);

fid = fopen("nozzle.fixa_outlet.tmp","w");
for k=1:length(outlet)-1
  node = outlet(k);
  fprintf(fid,"%d %d  %f\n",node,3,0.);
  fprintf(fid,"%d %d  %f\n",node,4,p_ref);
endfor
fclose(fid);

nnod = rows(xnod);
uini = [rho_ref u_ref v_ref p_ref];
uini = uini(ones(nnod,1),:);
asave("nozzle.ini.tmp",uini);
