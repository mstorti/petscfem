source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 hratio 1 1 hratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

[numnp,ndm]=size(xnod);
[numel,nen]=size(icone);

asave("poise2fases.nod.tmp",xnod);
asave("poise2fases.con.tmp",icone);

x=xnod(:,1);
y=xnod(:,2);
xmin = min(xnod);
xmax = max(xnod);

tol=1e-5;

nodes_in  = find(abs(x-xmin(1))<tol);
nodes_out = find(abs(x-xmax(1))<tol);
nodes_bot = find(abs(y-xmin(2))<tol);
nodes_top = find(abs(y-xmax(2))<tol);

% fijamos las velocidades de ambas fases
fixa_l = [];
nodes = nodes_in;
compo = 3*ones(length(nodes),1);
value = 1*ones(length(nodes),1);
fixa_l = [ fixa_l ; [nodes compo value ]];

nodes = nodes_in;
compo = 4*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_l = [ fixa_l ; [nodes compo value ]];
 
nodes = nodes_bot;
compo = 3*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_l = [ fixa_l ; [nodes compo value ]];

nodes = nodes_bot;
compo = 4*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_l = [ fixa_l ; [nodes compo value ]];

nodes = nodes_top;
compo = 3*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_l = [ fixa_l ; [nodes compo value ]];

nodes = nodes_top;
compo = 4*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_l = [ fixa_l ; [nodes compo value ]];

fixa_g = [];
nodes = nodes_in;
compo = 5*ones(length(nodes),1);
value = 1*ones(length(nodes),1);
fixa_g = [ fixa_g ; [nodes compo value ]];

nodes = nodes_in;
compo = 6*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_g = [ fixa_g ; [nodes compo value ]];
 
nodes = nodes_bot;
compo = 5*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_g = [ fixa_g ; [nodes compo value ]];

nodes = nodes_bot;
compo = 6*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_g = [ fixa_g ; [nodes compo value ]];

nodes = nodes_top;
compo = 5*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_g = [ fixa_g ; [nodes compo value ]];

nodes = nodes_top;
compo = 6*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa_g = [ fixa_g ; [nodes compo value ]];

fixa = [fixa_l ; fixa_g];

% fijamos la presion a la entrada y a la salida
if 0,
nodes = nodes_in;
compo = 1*ones(length(nodes),1);
value = 1*ones(length(nodes),1);
fixa = [ fixa ; [nodes compo value ]];
endif

nodes = nodes_out;
if 0,
nn    = length(nodes);
nodes([1,nn])=[];
endif
compo = 1*ones(length(nodes),1);
value = p_out*ones(length(nodes),1);
fixa = [ fixa ; [nodes compo value ]];

% fijamos la velocidad transversal a la salida del gas y el liquido
nodes = nodes_out;
compo = 4*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa = [ fixa ; [nodes compo value ]];

compo = 6*ones(length(nodes),1);
value = 0*ones(length(nodes),1);
fixa = [ fixa ; [nodes compo value ]];

% fijamos la fraccion de gas a la entrada
nodes = nodes_in;
compo = 2*ones(length(nodes),1);
value = alpha_g*ones(length(nodes),1);
fixa = [ fixa ; [nodes compo value ]];

asave("poise2fases.fixa.tmp",fixa);
fid = fopen("poise2fases.fixa_all.tmp","w");
if 0,
for j=1:numnp;
#  fprintf(fid,"%d %d %f\n",j,2,alpha_g);
#  fprintf(fid,"%d %d %f\n",j,3,u_l);
#  fprintf(fid,"%d %d %f\n",j,4,v_l);
#  fprintf(fid,"%d %d %f\n",j,5,u_g);
#  fprintf(fid,"%d %d %f\n",j,6,v_g);
  fprintf(fid,"%d %d %f\n",j,7,k);
  fprintf(fid,"%d %d %f\n",j,8,eps);
endfor
else
  XX=aload(save_file);
  p_x      = XX(:,1);
  alpha_g_x = XX(:,2);
  u_l_x    = XX(:,3);
  v_l_x    = XX(:,4);
  u_g_x    = XX(:,5);
  v_g_x    = XX(:,6);

for j=1:numnp;

if fixa_fase==1,
  fprintf(fid,"%d %d %f\n",j,1,p_x(j));
  fprintf(fid,"%d %d %f\n",j,3,u_l_x(j));
  fprintf(fid,"%d %d %f\n",j,4,v_l_x(j));
elseif fixa_fase==2,
  fprintf(fid,"%d %d %f\n",j,2,alpha_g_x(j));
  fprintf(fid,"%d %d %f\n",j,5,u_g_x(j));
  fprintf(fid,"%d %d %f\n",j,6,v_g_x(j));
elseif fixa_fase==0,
  ## disp(' AMBAS FASES INTERACTUANDO ')
else
					 error (" elegi mal la fase ")					 
endif
  fprintf(fid,"%d %d %f\n",j,7,k);
  fprintf(fid,"%d %d %f\n",j,8,eps);
endfor
endif
fclose(fid);

uini = [0.0 alpha_g 0.e-5 0.e-6 0.e-5 0.e-6 0.1 0.1];
uini=uini(ones(rows(x),1),:);
#uini=0.01*randn(rows(x),8);
asave("poise2fases.ini.tmp",uini);

bcconv = [];
nodes=nodes_bot;
N=length(nodes);
[bas,ibas]=sort(xnod(nodes,1));
nodes=nodes(ibas);
bcconv = [bcconv;[nodes(1:N-1) nodes(2:N)]];

nodes=nodes_out;
N=length(nodes);
[bas,ibas]=sort(xnod(nodes,2));
nodes=nodes(ibas);
bcconv = [bcconv;[nodes(1:N-1) nodes(2:N)]];

nodes=nodes_top;
N=length(nodes);
[bas,ibas]=sort(-xnod(nodes,1));
nodes=nodes(ibas);
bcconv = [bcconv;[nodes(1:N-1) nodes(2:N)]];

nodes=nodes_in;
N=length(nodes);
[bas,ibas]=sort(-xnod(nodes,2));
nodes=nodes(ibas);
bcconv = [bcconv;[nodes(1:N-1) nodes(2:N)]];

asave("poise2fases.bcconv.tmp",bcconv);

