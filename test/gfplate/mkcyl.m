## $Id: mkcyl.m,v 1.4 2005/01/24 22:01:13 mstorti Exp $
source("data.m.tmp");

w = zhomo([log(R) log(Rext) 0 pi],Nr+1,Nphi+1);
w = exp(w);
## w = 0.5*(w + 1./w); ## for a plate
[xnod,icone] = pfcm2fem(w);
nnod = rows(xnod);
icone = icone(:,[1 4 3 2]);

## slip on axis upstream
nline = Nphi+1;
upstream = (Nphi+1)*(1:Nr)';
pffixa("cylabso.fixa-ups.tmp",upstream,3);

## slip on axis downsstream
downstream = 1+(Nphi+1)*(0:Nr-1)';
pffixa("cylabso.fixa-down.tmp",downstream,3);

## Slip on skin
fid = fopen("cylabso.fixa-skin.tmp","w");
for k=1:Nphi+1
  x = xnod(k,:);
  fprintf(fid,"%g %d %d    %g %d %d\n",
	  xnod(k,1),k,2, xnod(k,2),k,3);
endfor
fclose(fid);

## Absorbing elements on exterior boundary

## Add two layers more of nodes (fictitious)
xfic = xnod(nline*Nr+(1:nline),:);
xnod = [xnod;
	xfic;
	xfic];

abso_con = [];
lm_nodes = [];
rs_nodes = [];
fid = fopen("cylabso.abso-con.tmp","w");
fid2 = fopen("cylabso.fixa-ext-std.tmp","w");
for k=1:Nphi+1
  ## real nodes (exterior first)
  rnodes = nline*(Nr-(0:2))+k;

  ## Make mesh equispaced in the normal direction
  ## at the absorbing boundary
  n1 = rnodes(1);
  n2 = rnodes(2);
  n3 = rnodes(3);
  xnod(n2,:) = 0.5*(xnod(n1,:)+xnod(n3,:));

  ## First fictitious node (lag. mult.)
  lagmulnd = rnodes(1)+nline;
  lm_nodes = [lm_nodes;lagmulnd];
  ## Second fictitious node (reference state)
  refstnode = lagmulnd + nline;
  rs_nodes = [rs_nodes; refstnode];
  abso_con = [abso_con;
	      rnodes,lagmulnd,refstnode];
  node = rnodes(1);
  nor = xnod(node,:);
  nor = nor/l2(nor);
  fprintf(fid,"%d %d %d %d %d   %g %g\n",
	  rnodes,lagmulnd,refstnode,nor);
  if nor(1)<=0
    ## Incoming flow. Impose rho,u,v
    fprintf(fid2,"%d %d %g\n",node,1,rhoref);
    fprintf(fid2,"%d %d %g\n",node,2,uref);
    fprintf(fid2,"%d %d %g\n",node,3,0.);
  else
    fprintf(fid2,"%d %d %g\n",node,4,pref);
  endif
endfor
fclose(fid);
fclose(fid2);

asave("cylabso.nod.tmp",xnod);
asave("cylabso.con.tmp",icone);

Uref = [rhoref,uref,0,pref];
pffixa("cylabso.fixa-ref.tmp",rs_nodes,1:4,Uref);

pffixa("cylabso.fixa-lm-nodes.tmp",lm_nodes,1:4);

uini = Uref;
uini(2) = uini(2);
nnod2 = rows(xnod);
uini = uini(ones(nnod2,1),:);
uini(lm_nodes,:) = 0;
real_nodes = complement(lm_nodes,(1:nnod2));
real_nodes = complement(rs_nodes,real_nodes)';
## uini(real_nodes,2) + uini(real_nodes,2) + 0.05;
asave("cylabso.ini.tmp",uini);

some = create_set([1:Nphi+1,1:nline:nnod,nline-1+(1:nline:nnod)])';
asave("cylabso.some-nodes.tmp",some);
