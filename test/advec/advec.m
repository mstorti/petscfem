##__INSERT_LICENSE__
## $Id: advec.m,v 1.7 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],N+1,M+1);
[xnod,icone]=pfcm2fem(w);
if !exist("noise")
  noise = 0.;
endif
hx = Lx/N;
hy = Ly/M;
nnod = rows(xnod);
if noise
  dx = noise*(2*rand(size(xnod))-1)*diag([hx hy]);
  tol = 1e-2*min([hx hy]);
  boundary = (abs(xnod(:,1))<tol | abs(xnod(:,2))<tol |
	      abs(xnod(:,1)-Lx)<tol | abs(xnod(:,2)-Ly)<tol);
  dx = leftscal(!boundary,dx);
  xnod = xnod + dx;
endif

icone=[icone(:,[1 4 3]);
       icone(:,[3 2 1])];
	
asave("advec.nod.tmp",xnod);

if strcmp(case_name,"rotating_cone")
  xele = pfnd2ele(xnod,icone,xnod);
  uele = [-Omega*(xele(:,2)-Ly/2) +Omega*(xele(:,1)-Lx/2)];
  fid = fopen("advec.con.tmp","w");
  for k=1:rows(icone)
    fprintf(fid,"%d %d %d    %f %f\n",icone(k,:),uele(k,:));
  endfor
  fclose(fid);
else 
  asave("advec.con.tmp",icone);
endif

if cone

  fid = fopen("advec.fixa.tmp","w");
  if !strcmp(case_name,"rotating_cone")
    for k=1:M+1
      node = k;
      fprintf(fid,"%d %d %f\n",node,1,0);
    endfor
  else
    for k=1:M+1
      node = k;
      fprintf(fid,"%d %d %f\n",node,1,0);
      fprintf(fid,"%d %d %f\n",node+N,1,0);
    endfor
    for k=2:M
      node = (k-1)*(M+1)+1;
      fprintf(fid,"%d %d %f\n",node,1,0);
      node = (k-1)*(M+1)+M+1;
      fprintf(fid,"%d %d %f\n",node+M,1,0);
    endfor
  endif
  fclose(fid);
  
  rcone = 0.4; sigma=0.2; 
  r = sqrt((xnod(:,1)-xini).^2+(xnod(:,2)-yini).^2);
  phi = exp(-r.^2/sigma^2);
  phi_c = exp(-rcone.^2/sigma^2);# phi on the cone border
  phi = phi-phi_c;
  phi = phi.*(phi>0.);
  asave("advec.ini.tmp",phi);

  uy>=0 || error("not uy<0 allowed");
  fid = fopen("advec.fixa-y0.tmp","w");
  if uy>0
    for k=2:N+1
      node = (M+1)*(k-1)+1;
      fprintf(fid,"%d %d %f\n",node,1,0);
    endfor
  endif  
  fclose(fid);

else

  if !exist("ydisc"); ydisc = .5; endif

  yydisc = ydisc*ux/sqrt(ux^2+uy^2);

  fid = fopen("advec.fixa.tmp","w");
  for k=1:M+1
    node = k;
    ## Coordinate orthogonal to velocity
    yy = xnod(node,:)*[-uy ux]'/sqrt(ux^2+uy^2);
    phi = tanh((yy-yydisc)/delta);
    fprintf(fid,"%d %d %f\n",node,1,phi);
  endfor
  fclose(fid);

  uy>=0 || error("not uy<0 allowed");
  fid = fopen("advec.fixa-y0.tmp","w");
  if uy>0
    for k=2:N+1
      node = (M+1)*(k-1)+1;
      ## Coordinate orthogonal to velocity
      yy = xnod(node,:)*[-uy ux]'/sqrt(ux^2+uy^2);
      phi = tanh((yy-yydisc)/delta);
      fprintf(fid,"%d %d %f\n",node,1,phi);
    endfor
  endif  
  fclose(fid);

endif
