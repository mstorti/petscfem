source("data.m.tmp");

w=zhomo([0 L 0 L],Nx+1,Ny+1,[1 0 1 1 0 1]);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);
nnod=rows(xnod);

zb=zeros(nnod,1);
xnod=[xnod(:,1:2) zb];
asave("plano.nod.tmp",xnod);
asave("plano.nod2d.tmp",xnod(:,1:2));
asave("plano.con.tmp",icone);
asave("plano.con0.tmp",icone(:,[1 4 2 3])-1);


if (0)
  if (turbulent)
    ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_ini*ones(nnod,1),\
	 kappa*ones(nnod,1),epsilon*ones(nnod,1)];
  else
    ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_ini*ones(nnod,1)];
  endif
endif

h_ini=base+ampli.*(1./(2*pi).*exp(-A*((xnod(:,1)-x0)./sigma).^2-B*((xnod(:,2)-y0)./sigma).^2));
ini=[u_ini*ones(nnod,1),v_ini*ones(nnod,1),h_ini];
asave("plano.ini.tmp",ini);

tol=1.e-4;
in=(1:Ny+1)';
out=((Nx)*(Ny+1)+1:(Nx+1)*(Ny+1))';
## l=(2*(Ny+1):Ny+1:(Nx)*(Ny+1))';
## r=((Ny+1)+1:Ny+1:(Nx)*(Ny+1))';
l=((Ny+1):Ny+1:(Nx+1)*(Ny+1))';
r=(1:Ny+1:(Nx+1)*(Ny+1))';

pffixa("plano.wall-fixa.tmp",[in;out;r;l],[1 2],[u_out,v_out]);
#pffixa("plano.wall-fixa.tmp",[in;out;r;l],[3],[h_out]);

if (turbulent)
  pffixa("plano.ke-fixa.tmp",(1:nnod)',[4 5],[kappa,epsilon]);
endif
if (0)
  pffixa("plano.in-fixa.tmp",in,[1 2],[u_in,v_in]);
  pffixa("plano.out-fixa.tmp",out,[3],[h_out]);
  pfconstr("plano.constr-slip-l.tmp",[l,l],1:2,[0. 1.]);
  pfconstr("plano.constr-slip-r.tmp",[r,r],1:2,[0. -1.]);
endif
