global Rint Rext Rext2 L Rmean

source("data.m.tmp");
rem(Ntheta,8)==0 || error("Ntheta must multiple of 8");
Rmean = (Rint+Rext)/2;

XNOD = [1 Rint*[cos(pi/4)  sin(pi/4)];
	2 Rmean*[cos(pi/4) sin(pi/4)];
	3 Rmean*cos(pi/4)  Rext;
	4 L                Rmean*sin(pi/4);
	5 L                Rext;
	6 Rint*[-cos(pi/4) sin(pi/4)];
	7 Rmean*[-cos(pi/4) sin(pi/4)];
	8 Rext*[-cos(pi/4) sin(pi/4)];
	9 Rint 0;
	10 Rmean 0;
	11 L 0;
	12 Rint*[cos(pi/4)  -sin(pi/4)];
	13 Rmean*[cos(pi/4) -sin(pi/4)];
	14 Rmean*cos(pi/4)  -Rext;
	15 L                -Rmean*sin(pi/4);
	16 L                -Rext;
	17 Rint*[-cos(pi/4) -sin(pi/4)];
	18 Rmean*[-cos(pi/4) -sin(pi/4)];
	19 Rext*[-cos(pi/4) -sin(pi/4)]	
	20 L  Rext2;
	21 Rmean*cos(pi/4)  Rext2;
	22 Rext2*[-cos(pi/4) sin(pi/4)];
	23 Rext2*[-cos(pi/4) -sin(pi/4)];
	24 Rmean*cos(pi/4)  -Rext2;
	25 L                -Rext2];

XNOD = XNOD(:,2:3);

ICONE = [1 2 7 6;
	 2 3 8 7;
	 2 4 5 3;
	 9 10 2 1;
	 10 11 4 2;
	 12 17 18 13;
	 13 18 19 14;
	 12 13 10 9;
	 13 15 11 10;
	 14 16 15 13;
	 6 7 18 17;
	 7 8 19 18];

ICONEE = [5 20 21 3;
	 3 21 22 8;
	 8 22 23 19;
	 19 23 24 14;
	 14 24 25 16];

H = [1 2 Nr/2;
     2 3 Nr/2;
     6 17 Ntheta/4;
     1 6 Ntheta/4;
     17 12 Ntheta/4;
     9 1 Ntheta/8;
     9 12 Ntheta/8;
     10 11 Nx;
     5 20 Next;
     20 21 Nx;
     21 22 Ntheta/4;
     22 23 Ntheta/4;
     23 24 Ntheta/4;
     24 25 Nx];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H);
external = mesher_bound(mesh,[5 3 8 19 14 16]);
outlet = mesher_bound(mesh,[16 15 11 4 5]);
skin = mesher_bound(mesh,[9 1 6 17 12 9]);

## Add a fictitious node for the constraints
xnod = [xnod;
	0 0];

asave("cylin.nod.tmp",xnod);
asave("cylin.con.tmp",icone);

nnod = rows(xnod);
uini = [1 0 0];
uini = uini(ones(nnod,1),:);
uini(:,1) = uini(:,1) + du_ini_pert * sin(pi/2*xnod(:,2)/Rext);
asave("cylin.ini.tmp",uini);

next = length(external);
normal = xnod(external(3:next),:) - xnod(external(1:next-2),:);
normal = [-normal(:,2) normal(:,1)];
normal = [0 1;
	  normal;
	  0 -1];
normal = leftscal(1./l2(normal),normal);

next = length(external);
normal = xnod(external(3:next),:) - xnod(external(1:next-2),:);
normal = [normal(:,2) -normal(:,1)];
normal = leftscal(1./l2(normal),normal);
normal = [0 1;
	  normal;
	  0 -1];

fid = fopen("cylin.normal.tmp","w");

for k=2:length(normal)-1
  node = external(k);
  fprintf(fid,"%f   %d %d      %f   %d %d    %f   %d %d\n",
	  normal(k,1),node,1, normal(k,2),node,2,
	  -normal(k,1),nnod,1);
endfor
fclose(fid);

fid = fopen("cylin.skin.tmp","w");
## Fix all values in the fictitious node
fprintf(fid,"%d %d    %f\n",nnod,1,1.);
fprintf(fid,"%d %d    %f\n",nnod,2,0.);
fprintf(fid,"%d %d    %f\n",nnod,3,0.);

for k=1:length(skin)
  node = skin(k);
  fprintf(fid,"%d %d    %f\n",node,1,0.);
  fprintf(fid,"%d %d    %f\n",node,2,0.);
endfor  
fclose(fid);

fid = fopen("cylin.outlet.tmp","w");
for k=1:length(outlet)
  node = outlet(k);
  fprintf(fid,"%d %d    %f\n",node,2,0.);
  fprintf(fid,"%d %d    %f\n",node,3,0.);
endfor  
fclose(fid);

[xnode,iconee,meshe] = mesher(XNOD,ICONEE,H);
external2 = mesher_bound(meshe,[5 3 8 19 14 16]);
external3 = mesher_bound(meshe,[20 21 22 23 24 25]);
out_up = mesher_bound(meshe,[5 20]);
out_down = mesher_bound(meshe,[25 16]);

fid = fopen("ext.fixa_ext.tmp","w");
for k=1:length(external2)
  node = external2(k);
  fprintf(fid,"%d %d   %f\n",node,1,1.);
endfor
fclose(fid);

asave("ext.nod.tmp",xnode);
asave("ext.con.tmp",iconee);

asave("ext.coupling_nodes.tmp",[external external2]);
