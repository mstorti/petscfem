### $Id: mkmesh.m,v 1.1 2002/10/09 13:38:44 mstorti Exp $ 

source("./data.m.tmp");

w=zhomo([R Rint 0 H],Nr+1,Nz+1,[1 0 1 yratio 0 1]);

[xnod,icone]=pfcm2fem(w);
[numnp,ndm]=size(xnod);
[numel,nen]=size(icone);

xmin = min(xnod);
xmax = max(xnod);
tol  = 1e-6;
nodes_in  = find(abs(xnod(:,1)-xmin(1))<tol);
nodes_out = find(abs(xnod(:,1)-xmax(1))<tol);
nodes_top = find(abs(xnod(:,2)-xmax(2))<tol);
nodes_bot = find(abs(xnod(:,2)-xmin(2))<tol);

if mesh_triangle==1
  icone=quad2tri_mesh(icone);
endif

[x3d,ico3d]= extrude (xnod,icone,nlay,dz);

if mesh_triangle==1
  ico3d = [ico3d(:,1:3),ico3d(:,3),ico3d(:,4:6),ico3d(:,6)];
endif

## input pressure
vaux = [];
nodes = (1:numnp)';
vaux = [nodes , (ndim+1)*ones(length(nodes),1) , pin*ones(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_in.tmp",vaux);

## output pressure
vaux = [];
nodes = numnp*(nlay)+(1:numnp)';
vaux = [nodes , (ndim+1)*ones(length(nodes),1) , pout*ones(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_out.tmp",vaux);

## non-slip velocity at top and bottom
if 1,
vaux = [];
nodes_base =  [nodes_in;nodes_top;nodes_bot;nodes_out];
for k=1:nlay+1
  nodes = nodes_base+(k-1)*numnp;
  vaux = [vaux ;[nodes , 1*ones(length(nodes),1) , 0*ones(length(nodes),1) \
	]];
  vaux = [vaux ;[nodes , 2*ones(length(nodes),1) , 0*ones(length(nodes),1) \
      ]]; 
  vaux = [vaux ;[nodes , 3*ones(length(nodes),1) , 0*ones(length(nodes),1) \
	]];
endfor

asave("bcconv_test_3d.fixa_wall.tmp",vaux);
endif

if weak_form==1,
## bcconv 
vaux=[];

if mesh_triangle==1,
## para triangulos anda 
icone_in = icone;
vaux = icone_in;
icone_out = icone(:,[3,2,1])+numnp*nlay;
vaux = [vaux; [ icone_out ] ];
else
## para quads  
icone_in = icone(:,[1,2,3,4]);
vaux = icone_in;
icone_out = icone(:,[4,3,2,1])+numnp*nlay;
vaux = [vaux; [ icone_out ] ];
endif

asave("bcconv_test_3d.bcconv.tmp",vaux);

endif

## write mesh
asave("bcconv_test_3d.nod.tmp",x3d);
asave("bcconv_test_3d.con.tmp",ico3d);

uini = ones(numnp*(nlay+1),1)*[0,0,0,pout];
asave("bcconv_test_3d.ini.tmp",uini);


return
## =================================================



## 3D extension for axisymmetric computation
angle = 1*pi/180;
rota = [cos(angle) , 0 , -sin(angle); 0,1,0; sin(angle), 0 , cos(angle)];
x3d = [xnod , xnod(:,1)*0];
x3d = [x3d; (rota*x3d')'];
icone = [icone , icone+numnp];

xmin = min(xnod);
xmax = max(xnod);
tol  = 1e-6;
nodes_in  = find(abs(xnod(:,2)-xmin(2))<tol);
nodes_out = find(abs(xnod(:,2)-xmax(2))<tol);
nodes_wall= find(abs(xnod(:,1)-xmax(1))<tol);
nodes_axis= find(abs(xnod(:,1)-xmin(1))<tol);

## input pressure
vaux = [];
nodes = nodes_in;
vaux = [nodes , (ndim+1)*ones(length(nodes),1) , pin*ones(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_in.tmp",vaux);
## output pressure
vaux = [];
nodes = nodes_out;
vaux = [nodes , (ndim+1)*ones(length(nodes),1) , pout*ones(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_out.tmp",vaux);
## non-slip wall velocity
vaux = [];
nodes = nodes_wall;
vaux = [nodes , 1*ones(length(nodes),1) , zeros(length(nodes),1) \
	];
vaux = [vaux ;[nodes , 2*ones(length(nodes),1) , zeros(length(nodes),1) \
	       ]]; 
asave("bcconv_test_3d.fixa_wall.tmp",vaux);
## slip axis velocity
vaux = [];
nodes = nodes_axis;
vaux = [nodes , 1*ones(length(nodes),1) , zeros(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_axis.tmp",vaux);

## w fixed for all nodes 
vaux = [];
nodes = (1:numnp)';
vaux = [nodes , 3*ones(length(nodes),1) , zeros(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_w_all.tmp",vaux);

if 0,
## output radial velocity
vaux = [];
nodes = nodes_out;
vaux = [nodes , 1*ones(length(nodes),1) , zeros(length(nodes),1) \
	];
asave("bcconv_test_3d.fixa_uout.tmp",vaux);
endif

## periodic constrain for axisymmetric computation
## Imposing axisymmetric fields 
uno   = ones(numnp,1);
node1 = numnp+(1:numnp)';
node2 = (1:numnp)';
no2   = [0;0;1];
no1   = rota*no2;

vaux = [];
vaux = [uno, node1 ,(ndim+1)*uno , -uno, node2, (ndim+1)*uno];
vaux = [vaux; [uno, node1 ,2*uno , -uno , node2 , 2*uno]];
asave("bcconv_test_3d.axisym_p_v.tmp",vaux);

vaux_vel = [ ...
	[uno, node1 ,1*uno , -uno*rota(1,1) , node2 , 1*uno, -uno*rota(1,3) , node2 , 3*uno];
	[uno, node1 ,3*uno , -uno*rota(3,1) , node2 , 1*uno, -uno*rota(3,3) , node2 , 3*uno]];
asave("bcconv_test_3d.axisym_u_w.tmp",vaux_vel);

## write mesh

asave("bcconv_test_3d.nod.tmp",x3d);
asave("bcconv_test_3d.con.tmp",icone);

uini = ones(2*numnp,1)*[uin,0,0,pin];
##uini(:,1:3)=0.01*rand(2*numnp,ndim);
asave("bcconv_test_3d.ini.tmp",uini);

## bcconv 
vaux=[];

nodes = nodes_in;
N = length(nodes);
idir = (1:N-1);
tmp1 = [nodes(idir),nodes(idir+1)];
vaux = [vaux;[tmp1 numnp+tmp1(:,[2,1])]];

nodes = nodes_out;
N = length(nodes);
idir = (1:N-1);
tmp2 = [nodes(idir+1),nodes(idir)];
vaux = [vaux;[tmp2 numnp+tmp2(:,[2,1])]];

asave("bcconv_test_3d.bcconv.tmp",vaux);

