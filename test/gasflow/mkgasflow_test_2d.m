##__INSERT_LICENSE__
## $Id: mkgasflow_test_2d.m,v 1.3 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");

#w=zhomo([0 D 0 H],Nr+1,Nz+1,[1 rratio 1 1 hratio 1.2*hratio]);
w=zhomo([0 D 0 H],Nr+1,Nz+1);

[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

[numnp,ndm]=size(xnod);
[numel,nen]=size(icone);

asave("gasflow_test_2d.nod.tmp",xnod);
asave("gasflow_test_2d.con.tmp",icone);

x=xnod(:,1);
y=xnod(:,2);
xmin = min(xnod);
xmax = max(xnod);

tol=1e-5;

nodes_lef = find(abs(x-xmin(1))<tol);
nodes_rig = find(abs(x-xmax(1))<tol);
nodes_bot = find(abs(y-xmin(2))<tol);
nodes_top = find(abs(y-xmax(2))<tol);

## Fixations
## inlet
fixa = [];
nodes = nodes_lef;
values = [rho_ref,u_ref,v_ref,p_ref];
for icompo=[1,2,3]
  compo = icompo*ones(length(nodes),1);
  value = values(icompo)*ones(length(nodes),1);
  fixa  = [ fixa ; [nodes compo value ]];
endfor
asave("gasflow_test_2d.fixa_in.tmp",fixa);

## outlet
fixa = [];
nodes = nodes_rig;
values = [rho_ref,u_ref,0,p_ref];
for icompo=[4]
  compo = icompo*ones(length(nodes),1);
  value = values(icompo)*ones(length(nodes),1);
  fixa  = [ fixa ; [nodes compo value ]];
endfor
asave("gasflow_test_2d.fixa_out.tmp",fixa);

## top and bottom
fixa = [];
nodes = [nodes_top;nodes_bot];
values = [rho_ref,u_ref,v_ref,p_ref];
for icompo=[3]
  compo = icompo*ones(length(nodes),1);
  value = values(icompo)*ones(length(nodes),1);
  fixa  = [ fixa ; [nodes compo value ]];
endfor
asave("gasflow_test_2d.fixa_wall.tmp",fixa);

## Bcconv
if weak_form,
bcconv = [];
nodes=nodes_lef;
N=length(nodes);
[bas,ibas]=sort(xnod(nodes,1));
nodes=nodes(ibas);
bcconv = [bcconv;[ nodes(2:N) nodes(1:N-1)]];
nodes=nodes_rig;
N=length(nodes);
[bas,ibas]=sort(xnod(nodes,1));
nodes=nodes(ibas);
bcconv = [bcconv;[ nodes(1:N-1) nodes(2:N) ]];
asave("gasflow_test_2d.bcconv.tmp",bcconv);
endif

## Initial condition
uini = [rho_ini,u_ini,v_ini,p_ini];
uini=uini(ones(rows(x),1),:);
# pressure perturbation
DP     = 0.001*p_ini;
sigma  = 0.05;
xc     = (max(xnod)+min(xnod))/2;
r      = sqrt(sum((xnod-ones(numnp,1)*xc)'.^2))';
uini(:,4) = uini(:,4)+DP*exp(-r.^2/sigma);
asave("gasflow_test_2d.ini.tmp",uini);
