##__INSERT_LICENSE__
## $Id: les.m,v 1.5 2003/03/13 17:03:27 mstorti Exp $
source("data.m.tmp");
nlay=nz;			# number of element layers
Lz=1;				# length in the z direction
                                # Side Lx is 1.
alpha=1.7;			# Controls refinement towards solid
				# wall. alpha=1: linear, =2: quadratic

#########

## basic 2d mesh
w=zhomo([0 1 0 1],nx+1,nx+1);
[xnod,icone]=pfcm2fem(w);
clear w
icone=icone(:,[1 4 3 2]);

dz=Lz/2/nlay;
nnol=(nx+1)^2; # number of nodes in a layer
ic31=[icone icone+nnol]; # con. table of a slab layer
x3=[];
ic3=[];

## On the top (z=Lz)
fixa = [(1:nnol)'+nnol*nlay 1*ones(nnol,1) ones(nnol,1); # u fixed to 1
        (1:nnol)'+nnol*nlay 2*ones(nnol,1) zeros(nnol,1); # v,w fixed to 0.
        (1:nnol)'+nnol*nlay 3*ones(nnol,1) zeros(nnol,1)]; 

## Side x=1 periodic with x=0
perix = (1:nx+1)';
perix=[perix perix+(nx+1)*nx];

## Side y=1 periodic with y=0
periy = (1:nx+1:(nx+1)^2)';
periy=[periy periy+nx];

perix3=[];
periy3=[];

## Constructs 3D con. table
for k=1:nlay
  ic3=[ic3;
       ic31+nnol*(k-1)];
endfor

## Constructs 3D node array and perio. table
for k=1:nlay+1
  x3=[x3;
      xnod,(k-1)*dz*ones(nnol,1)];
  perix3=[perix3;
          perix+nnol*(k-1)];
  periy3=[periy3;
          periy+nnol*(k-1)];
endfor

## refines towards solid wall
z=x3(:,3);
Lz2=Lz/2;
z=(z/Lz2).^alpha*Lz2;
x3(:,3) = z;
clear z

peri3=[perix3;
       periy3];
O=ones(rows(peri3),1);

## Adds coefficients for the constraints
peri=[-O peri3(:,1) O O peri3(:,2) O;
      -O peri3(:,1) 2*O O peri3(:,2) 2*O;
      -O peri3(:,1) 3*O O peri3(:,2) 3*O;
      -O peri3(:,1) 4*O O peri3(:,2) 4*O];

## Initial state vector
nnod = rows(x3);
z=x3(:,3);
f=z.*(Lz-z)/Lz2^2;		# From 0 on the wall to one on the
				# symmetry plane (z=Lz/2)
x = x3(:,1);
y = x3(:,2);
z = x3(:,3);
u=[f+du_pert*sin(L.*sin(2*pi*y)) 0*z 0*z 0*z];
u(:,3) = u(:,3).*f;

## Fixations for the wall (if don't want to use LES)
fixaw = [(1:nnol)' ones(nnol,1) zeros(nnol,1); # set x comp.
         (1:nnol)' 2*ones(nnol,1) zeros(nnol,1); # set y comp. to 0.
         (1:nnol)' 3*ones(nnol,1) zeros(nnol,1)]; # set z comp. to 0.

## Save data to file
asave("les.nod.tmp",x3);
asave("les.con.tmp",ic3);
asave("les.fixa.tmp",fixa);
asave("les.fixaw.tmp",fixaw);
asave("les.peri.tmp",peri);
asave("les.wall.tmp",icone);
asave("les.inirand.tmp",u);
