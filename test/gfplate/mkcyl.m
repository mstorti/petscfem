## $Id: mkcyl.m,v 1.1 2005/01/24 00:48:32 mstorti Exp $
source("data.m.tmp");

w = zhomo([log(R) log(Rext) 0 pi],Nr+1,Nphi+1);
w = exp(w);
[xnod,icone] = pfcm2fem(w);
nnod = rows(xnod);
icone = icone(:,[1 4 3 2]);

## slip on axis upstream
nline = Nphi+1;
upstream = (Nphi+1)*(1:Nr+1)';
pffixa("cylabso.fixa-ups.tmp",upstream,3);

## slip on axis downsstream
downstream = 1+(Nphi+1)*(0:Nr)';
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
for k=1:Nphi+1
  ## real nodes (exterior first)
  rnodes = nline*(Nr-(0:2))+k;
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
endfor
fclose(fid);
## asave("cylabso.abso-con.tmp",abso_con);

asave("cylabso.nod.tmp",xnod);
asave("cylabso.con.tmp",icone);

Uref = [rhoref,uref,0,pref];
pffixa("cylabso.fixa-ref.tmp",rs_nodes,1:4,Uref);

uini = Uref;
nnod2 = rows(xnod);
uini = uini(ones(nnod2,1),:);
uini(lm_nodes,:) = 0;
asave("cylabso.ini.tmp",uini);
