source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

## add nodes for the stream
xnod = [xnod;
	xnod(1:Ny+1,:)];

if var_eta
  xe = pfnd2ele(xnod,icone,xnod(:,1));
  icone = [icone eta0+(etaL-eta0)/Lx*xe];
endif

## number of elements in the stream
nelst = Ny;

## nodes on the aquifer, under the stream
naq = (1:Ny+1)';

## nodes on the stream
nst = (Nx+1)*(Ny+1)+(1:Ny+1)';

## connectivities for the stream
icost = [nst(1:nelst) nst(2:nelst+1)];

## connectivities for the stream-loss elemset
icostl = [icost naq(1:nelst) naq(2:nelst+1)];

x=xnod(:,1);
y=xnod(:,2);

## arc length coordinate on the stream
s = xnod(:,2);

## Set bottom height for the channel
h_b = eta0+(etaL-eta0)/Lx*xnod(:,2);
h_b(nst) = -slope * s(nst);
xnod = [xnod h_b];

## Initial state
uini=zeros(rows(xnod),1);
uini(nst,1) = 0.3*ones(nelst+1,1);

asave("aquist.nod.tmp",xnod);
asave("aquist.con.tmp",icone);
asave("aquist.con.stream.tmp",icost);
asave("aquist.con.stream_loss.tmp",icostl);
asave("aquist.ini.tmp",uini);

tol=1e-5;

fid = fopen("aquist.fixa.tmp","w");
x0 = (1:Ny+1)';
for k=x0'
  fprintf(fid,"%d 1 0.\n",k);
endfor

xL = Nx*(Ny+1)+(1:Ny+1)';
for k=xL'
  fprintf(fid,"%d 1 1.\n",k);
endfor
fclose(fid);
