source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],N+1,M+1);
[xnod,icone]=pfcm2fem(w);

icone=[icone(:,[1 4 3]);
       icone(:,[3 2 1])];
	     
asave("advec.nod.tmp",xnod);
asave("advec.con.tmp",icone);

if strcmp(case_name,"cone") || strcmp(case_name,"rotating_cone")

  fid = fopen("advec.fixa.tmp","w");
  for k=1:M+1
    node = k;
    fprintf(fid,"%d %d %f\n",node,1,0);
  endfor
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
