##__INSERT_LICENSE__
## $Id: mksqcav.m,v 1.1 2005/01/07 01:44:20 mstorti Exp $
source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

w=zhomo([0 1 0 1],N+1,N+1,[1 hratio 1 1 hratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

if use_triangles
  icone = [icone(:,[1 2 3]);
	   icone(:,[3 4 1])];
endif

asave("sqcav.nod.tmp",xnod);
asave("sqcav.con.tmp",icone);

nnod = rows(xnod);
uini = [u_ini 0 0];
uini = uini(ones(nnod,1),:);
asave("sqcav.ini.tmp",uini);

x=xnod(:,1);
y=xnod(:,2);

tol=1e-5;
lid = find(abs(y-1)<tol)';
bot = find(abs(y)<tol)';

lat = find(abs(x-0)<tol | abs(x-1)<tol)';

## constant pressure on top
pffixa("sqcav.lid.tmp",complement(lat,lid),4);

## null velocity at all nodes except contact line
tmp = union(lat,bot);
tmp = complement(lid,tmp);
pffixa("sqcav.ns-wall.tmp",tmp,[1,2]);

## only horizontal velocity imposed at water-line
tmp = intersection(lat,lit);
pffixa("sqcav.ns-wl.tmp",tmp,1);

## null displacement at bottom
pffixa("sqcav.mm-bot.tmp",bot,[1,2]);

## null horizontal displacements at lateral
tmp = complement(bot,lat);
pffixa("sqcav.mm-bot.tmp",tmp,1);

## displacements imposed at free surface
tmp = complement(bot,lat);
pffixa("sqcav.mm-bot.tmp",tmp,1);
