## __INSERT_LICENSE__
## $Id: mkcloscyl.m,v 1.1 2002/08/17 22:27:47 mstorti Exp $

## Do not save the history commands when running in batch
if !index("program_invocation_name","/octave")
  saving_history = 0;
endif

source("./data.m.tmp");

rem(n,2)==0 || error("n should be pair, n= %d\n",n);

## build a square
w=zhomo([-1 1 -1 1],n+1,n+1);
[x2,ic2]=pfcm2fem(w);
nnod=rows(x2);
## add coordinate z=1
x2=[x2 ones(nnod,1)];

## paste two faces along x=1 and y=1
xnod=x2;
icone=ic2;
[xnod,icone]=paste(xnod,icone,x2(:,[2 3 1]),ic2);
[xnod,icone]=paste(xnod,icone,x2(:,[3 1 2]),ic2);

## paste the symmetry (3 more faces)
[xnod,icone]=paste(xnod,icone,-xnod,icone(:,[1 4 3 2]));
zelem = pfnd2ele(xnod,icone,xnod(:,3));

## find those elements with z>0
up = find(zelem>0);
icone = icone(up,:);

## purge nodes
[xnod,icone] = purgend(xnod,icone);

## Project on unit sphere
R=sqrt(sum((xnod.^2)'))';
xnod = leftscal(1./R,xnod);
clear R

## take spherical (Theta and phi) coordinates
rho = sqrt(sum((xnod(:,1:2).^2)'))';
Theta = atan2(rho,xnod(:,3));
phi = atan2(xnod(:,2),xnod(:,1));

## 2d surface
r = Theta/max(Theta);		# normalize to one
r = onedstr([rratio 0 1],r);	# refine toward walls
x2 = [r.*cos(phi) r.*sin(phi)];	# convert to cartesian coords.
xnod=x2;
clear x2

## scale to radius=RR
r = sqrt(sum((xnod.^2)'))';
xnod = xnod*RR;

## find nodes in the intersectin of free-surface and container
r = sqrt(sum((xnod.^2)'))';
waterline = find(abs(r-RR)<1e-4);

## order nodes so as to make a surface grid
Theta = atan2(xnod(:,2),xnod(:,1));
[bid,indx] = sort(Theta(waterline));
waterline =waterline(indx);
nplane = length(xnod);          # number of nodes in a plane
clear Theta

## Find node at the center 
center = find(sum(abs(xnod'))'<1e-6);
length(center)==1 || error("found more than one node at the center");

## extrude mesh in z direction
[x3d,i3d] = extrude(xnod,icone,nz,1/nz);
if use_tetra
  error("Should use the hexasplit.bin program");
  ##   i3d = hexasplit(i3d);
endif

## refine in z direction
z=x3d(:,3);
if exist("refine_to_bottom") && refine_to_bottom
  z=-Lz*onedstr([1 zratio 1],z);
else
  z=-Lz*onedstr([1 0 zratio],z);
endif
x3d=[x3d(:,1:2) z];

nnod = rows(x3d);		# number of real nodes
nodf_u = nnod+1;                # fictitious node for ux,uy
nodf_r = nnod+2;                # fictitious node for rotations
nfic = 2;			# number of fictitious nodes
x3d = [x3d;
       0 0 0;
       0 0 0];
nelem = rows(i3d);
if use_rot
  i3d = [i3d, nodf_u*ones(nelem,1), nodf_r*ones(nelem,1)];
endif

asave("cylinder.nod.tmp",x3d);
asave("cylinder.con.tmp",i3d);

## Connectivities for the surface elements for integration of free
## surface level
asave("cylinder.integr_top.tmp",icone + nnod);

## Connectivities for the bottom elements for integration of forces
asave("cylinder.bottom_con.tmp",[icone(:,[1 4 3 2])+nplane*nz ones(rows(icone),layers*4)]);
asave("cylinder.force_int_top.tmp",[icone ones(rows(icone),layers*4)]);

## list of fictitious nodes for dummy elemset
asave("cylinder.fic.tmp",[nodf_u; nodf_r]);

## nodes on the lateral wall
wall=[];
for k=0:nz;
  wall = [wall;
          k*nplane+waterline];
endfor
bottom = nz*nplane+(1:nplane)';
top = (1:nplane)';

## Connectivities of panels on the lateral wall
## (for force computation and bcconv)
fidbcc = fopen("cylinder.bcconv.tmp","w");
fid = fopen("cylinder.wall_panel_con.tmp","w");
nwl = length(waterline);
for j=1:nwl
  n1 = waterline(j);
  j1 = j+1;
  if j1>nwl, j1=1; endif
  n2 = waterline(j1);
  for k=1:nz
    conn = [((k-1)*nplane)+[n2 n1] (k*nplane)+[n1 n2] ];
    fprintf(fid,"%d %d %d %d",conn([1 4 3 2]));
    for j=1:layers, fprintf(fid," 1 1 1 1 "); endfor;
    fprintf(fid,"\n");
    fprintf(fidbcc,"%d %d %d %d\n",conn([1 4 3 2]));
  endfor
endfor
fclose(fid);
fclose(fidbcc);


## Impose slip condition on lateral wall
fid = fopen("cylinder.no_slip.tmp","w");
all_no_slip = create_set([wall;
			  bottom;
			  top]);
for k=1:length(all_no_slip)
  ## Impone velocidad igual a una cierta velocidad de rotacion
  node = all_no_slip(k);
  fprintf(fid,"%d 1 %g\n",all_no_slip(k),-Omega*x3d(node,2));
  fprintf(fid,"%d 2 %g\n",all_no_slip(k),+Omega*x3d(node,1));
  fprintf(fid,"%d 3 0.\n",all_no_slip(k));
endfor

## Impose to a given fixed value k and eps
tke = 0.1;
fidke = fopen("cylinder.ke_fix.tmp","w");
for k=1:nnod
  fprintf(fidke,"%d %d %f\n",k,ndim+2,tke);
  fprintf(fidke,"%d %d %f\n",k,ndim+3,tke);
endfor
fclose(fidke);

nnod = rows(x3d);
uini = [zeros(1,ndim+1), tke ,tke];
uini = uini(ones(nnod,1),:);
uini(:,1) = -Omega*x3d(:,2);
uini(:,2) = +Omega*x3d(:,1);
asave("cylinder.ini.tmp",uini);

## Pressure node
node_press = center+nz*nplane;
fprintf(fid,"%d 4 0.\n",node_press);
fclose(fid);

