##__INSERT_LICENSE__
## $Id: proc.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");
NN=(N+1);
NNy = Ny+1;
nnod = NN*NNy;
ndof=3;
xx=aload("wallke.nod.tmp");
x=xx(1:NNy:nnod,1);
if !exist("state") 
  args=[];
else
  args=split(state,"/");
  eval(["clear " state]);
endif

nargs = rows(args);
if nargs<3
  indx=1;
else
  indx = str2num (args(3,:));
endif

if nargs<2
  filename = "save.state";
else 
  filename = deblank(args(2,:));
  if length(filename)==0
    filename="outvector0.out";
  endif
endif

if nargs<1;
  ss="";
else
  ss=deblank(args(1,:));
endif

u=read_state(filename,nnod,ndof,indx);
out=NNy:NNy:nnod;
vout=u(out,2);
cntr=(NN-1)*NNy+(1:NNy)';
yc=xx(cntr,2);
uc=u(cntr,2);

U=reshape(u(:,1),NNy,NN)';
V=reshape(u(:,2),NNy,NN)';
P=reshape(u(:,3),NNy,NN)';

Q=sum(leftscal(diff(x),xcent(V)));

if exist("ss") && length(ss)>0
  eval([ss ".U = U;"]);
  eval([ss ".V = V;"]);
  eval([ss ".E = E;"]);
  eval([ss ".Q = Q;"]);
  clear ss
endif

