###key parti.m

## Reads a node partitioning and the connectivity and
## computes an element partitioning for visualization with DX

if 1
  npart=aload("node-part.dat");
  npart=npart(:,2);
  icone=aload("cubcav.con-tetra.tmp");
  xnod = aload("cubcav.nod.tmp");
endif

icosub=icone;
nelem=rows(icosub);
for k=1:columns(icosub)
  icosub(:,k) = npart(icosub(:,k));
endfor

nsubdo = max(npart)+1;
subdo = -ones(nelem,1);
vert = zeros(nelem,1);
for k=0:nsubdo-1
  vert1 = sum(icosub'==k)';
  subdo = choose(vert<vert1,k*ones(nelem,1),subdo);
endfor

for k=0:nsubdo-1
  printf("in subd %d: %d\n",k,sum(subdo==k));
endfor

nnod = rows(npart);
in_subdo = zeros(nnod,nsubdo);
asave("epart.dat",subdo);

xele=pfnd2ele(xnod,icone,xnod);

colors = [1 0 0;
	  0 1 0;
	  1 1 0;
	  0 1 1;
	  1 0 1];
zsubdo=zeros(nsubdo,1);

for k=0:nsubdo-1
  indx = find(subdo==k);
  zsubdo(k+1) = mean(xele(indx,3));
endfor

[bid,sbdindx] = sort(-zsubdo);

for k=1:nsubdo
  indx = find(subdo==sbdindx(k)-1);
  color=colors(rem(k,rows(colors))+1,:);
  color=color(ones(rows(indx),1),:);
  eval(["asave(\"icone" int2str(k-1) ".dat\",icone(indx,:)-1);"]);
  eval(["asave(\"colors" int2str(k-1) ".dat\",color);"]);
endfor
