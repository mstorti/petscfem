source("data.m");

## Physical parameters
uu=[u v];
Ax=u*eye(ndof);
Ay=v*eye(ndof);
DD=[D 0; 0 D];
Dxx=DD(1,1)*eye(ndof);
Dxy=DD(1,2)*eye(ndof);
Dyy=DD(2,2)*eye(ndof);
S=[s s s]';
RR=R*eye(ndof);
CP = Cp*eye(ndof);

if full_jacs==0
  ## this means no full_jacs
elseif full_jacs==1
  fid=fopen("fulljacdef.tmp","w");
  if strcmp(dif_type,"scalar_per_field")
    RR=diag(RR);
    RR=RR+R*Rfluc*(2*rand(size(RR))-1);
    fprintf(fid,"reactive_jacobians_type \"scalar_per_field\"\n");
    fprintf(fid,"reactive_jacobians");
    for j=1:ndof
      fprintf(fid," %f ",RR(j));
    endfor
    fprintf(fid,"\n\n");
    RR=diag(RR);
  endif

  if strcmp(dif_type,"scalar_per_field")
    fprintf(fid,"diffusive_jacobians_type \"scalar_per_field\"\n");
    fprintf(fid,"diffusive_jacobians");
    Dxx=0; Dxy=0; Dyy=0;
    for k=1:ndof
      ddd=D*(1+Dfluc*(2*rand-1));
      Dxx(k,k)=ddd;
      fprintf(fid," %f ",ddd);
    endfor
    fprintf(fid,"\n\n");
    Dyy=Dxx;
  endif
  
  if strcmp(adv_type,"vector_per_field")
    fprintf(fid,"advective_jacobians_type \"vector_per_field\"\n");
    fprintf(fid,"advective_jacobians");
    Ax=0; Ay=0;
    for k=1:ndof
      uuu=u*(1+ufluc*(2*rand-1));
      Ax(k,k)=uuu;
      vvv=v*(1+ufluc*(2*rand-1));
      Ay(k,k)=vvv;
      fprintf(fid," %f %f ",uuu,vvv);
    endfor
    fprintf(fid,"\n\n");
    Dyy=Dxx;
  endif
  
  if strcmp(enthalpy_type,"scalar_per_field")
    fprintf(fid,"enthalpy_jacobians_type \"scalar_per_field\"\n");
    fprintf(fid,"enthalpy_jacobians_type");
    CP=diag(Cp+Cpfluc*(2*rand(3,1)-1));
    for k=1:ndof
      fprintf(fid," %f ",CP(k,k));
    endfor
    fprintf(fid,"\n\n");
  endif
  
  if strcmp(source_type,"full")
                                #    fprintf(fid,"source_term_type \"full\"\n");
                                #    fprintf(fid,"source_term");
    S=0;
    for k=1:ndof
      sss=s*(1+sfluc*(2*rand-1));
      S(k)=sss;
                                #      fprintf(fid," %f ",sss);
    endfor
                                #    fprintf(fid,"\n\n");
  endif
  
  fclose(fid);
elseif full_jacs==2
  fid=fopen("fulljacdef.tmp","w");
  ## We should check that R is positive definite
  RR=R*eye(ndof)+R*Rfluc*(2*rand(size(RR))-1);
  fprintf(fid,"reactive_jacobians_type \"full\"\n");
  fprintf(fid,"reactive_jacobians");
  for j=1:ndof
    for k=1:ndof
      fprintf(fid," %f ",RR(j,k));
    endfor
  endfor
  fprintf(fid,"\n\n");

  CP=log(Cp)*(eye(ndof)+log((Cp-Cpfluc)/Cp)*(2*rand(size(CP))-1));
  CP=(CP+CP')/2;
  CP=expm(CP);
	      
  fprintf(fid,"enthalpy_jacobians_type \"full\"\n");
  fprintf(fid,"enthalpy_jacobians");
  for j=1:ndof
    for k=1:ndof
      fprintf(fid," %f ",CP(j,k));
    endfor
  endfor
  fprintf(fid,"\n\n");

  fprintf(fid,"diffusive_jacobians_type \"full\"\n");
  fprintf(fid,"diffusive_jacobians ");
  ## We should check that D is a positive
  ## definite operator
  Dxx=D*(eye(ndof)+Dfluc*(2*rand(ndof)-1));
  Dxx=(Dxx+Dxx')/2;
  Dxy=D*(          Dfluc*(2*rand(ndof)-1));
  Dxy=(Dxy+Dxy')/2;
  Dyy=D*(eye(ndof)+Dfluc*(2*rand(ndof)-1));
  Dyy=(Dyy+Dyy')/2;
  for j=1:ndof
    for k=1:ndof
      fprintf(fid,"\\\n   %f  %f  %f  %f   ", \
              Dxx(j,k),Dxy(j,k),Dxy(j,k),Dyy(j,k));
    endfor
  endfor
  fprintf(fid,"\n\n");
  
  fprintf(fid,"advective_jacobians_type \"full\"\n");
  fprintf(fid,"advective_jacobians ");
  Ax=u*(2*rand(ndof)-1);
  Ax=u*(eye(ndof)+ufluc*(Ax+Ax')/2);
  Ay=v*(2*rand(ndof)-1);
  Ay=v*(eye(ndof)+ufluc*(Ay+Ay')/2);
  for k=1:ndof
    fprintf(fid,"\\\n    %f %f %f  ",Ax(k,:));
  endfor
  for k=1:ndof
    fprintf(fid,"\\\n    %f %f %f  ",Ay(k,:));
  endfor
  
  if strcmp(source_type,"full")
                                #    fprintf(fid,"source_term_type \"full\"\n");
                                #    fprintf(fid,"source_term");
    S=0;
    for k=1:ndof
      sss=s*(1+sfluc*(2*rand-1));
      S(k)=sss;
                                #      fprintf(fid," %f ",sss);
    endfor
                                #    fprintf(fid,"\n\n");
  endif
  
  fclose(fid);
else
  full_jacs
  error("Value of full_jacs unknown");
endif 

## Number of nodes/elements in x,y
Nx=nx+1;
Ny=ny+1;

w=zhomo([0 Lx 0 Ly],nx+1,ny+1);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);

yy=xnod(1:ny+1,2);


## left to right peri
peri=[(ny+1)*nx+(1:ny+1)' (1:ny+1)'];
## upper to lower
nnod=(nx+1)*(ny+1);
peri=[peri;
      ny+(1:(ny+1):nnod)'  (1:(ny+1):nnod)'];
fid=fopen("newff.peri.tmp","w");
for j=1:ndof
  for k=1:rows(peri);
    fprintf(fid,"%d   %d   %d   %d   %d   %d\n", \
            -1,peri(k,1),j,1,peri(k,2),j);
  endfor
endfor
fclose(fid);

if use_bcconv
  ## bcconv
  nodes_right=(ny+1)*nx+(1:ny+1)';
  nodes_left=(ny+1:-1:1)';
  nodes_top=flipud(ny+(1:(ny+1):nnod)');
  nodes_bot=(1:(ny+1):nnod)';
#    bcc_elems=[nodl2elem(nodes_left);
#               nodl2elem(nodes_bot);
#               nodl2elem(nodes_right);
#               nodl2elem(nodes_top)];
  bcc_elems=nodl2elem(nodes_right);
  asave("newff.bcconv.tmp",bcc_elems);

  fid=fopen("newff.rfixa.tmp","w");
  for j=1:rows(nodes_left);
    for k=1:ndof;
      fprintf(fid,"   %d  %d  %f  \n", \
              nodes_left(j),k,rval);
    endfor
  endfor
  fclose(fid);

endif                           # { use_bcconv }


fixa=[];
for j=1:3
  fixa=[fixa;
        (1:ny+1)' j*ones(ny+1,1) ones(ny+1,1)];
endfor

fid=fopen("newff.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);

nele=rows(icone);
if per_elem_prop
  xe=zeros(nele,2);
  for k=1:4
    xe=xe+xnod(icone(:,k),:);
  endfor
  xe=xe/4;
endif

##sour=cos(kx*xe(:,1)*2*pi/Lx).*cos(ky*xe(:,2)*2*pi/Ly);
kwave=[kx*2*pi/Lx ky*2*pi/Ly];
phase=kwave(1)*xe(:,1)+kwave(2)*xe(:,2);
phase=exp(i*phase);
sour=real(phase)*S.';

asave("newff.nod.tmp",xnod);
                                #asave("newff.con.tmp",icone);
fid=fopen("newff.con.tmp","w");
for k=1:nele
  fprintf(fid,"%d %d %d %d",icone(k,1),icone(k,2),icone(k,3),icone(k,4));
  if per_elem_prop && !use_bcconv
    fprintf(fid," %f %f %f ",sour(k,1:3));
  endif
  fprintf(fid,"\n");
endfor
fclose(fid);

## bcconv en todo el fondo y la tapa
#  bcconv = [nx*Ny+(1:Ny)';
#            (nx*Ny:-Ny:Ny)'];

#  some = Nx*Ny;
#  asave("newff.some.tmp",some);

#  nbc=rows(bcconv);
#  bcconv = [bcconv(1:nbc-1) bcconv(2:nbc)];
#  asave("newff.bcconv.tmp",bcconv);

##--<*>---//---<*>---//---<*>---//---<*>---//---<*>---//
##--<*>---//--- ANALYTICAL SOLUTION -<*>---//---<*>---//
##--<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
kwu=kwave(1)*Ax+kwave(2)*Ay;
kDk = kwave(1)*kwave(1)*Dxx+2*kwave(1)*kwave(2)*Dxy+kwave(2)*kwave(2)*Dyy;
beta = CP \ (RR + i*kwu + kDk);
phase=kwave(1)*xnod(:,1)+kwave(2)*xnod(:,2);
phase=exp(i*phase);
if !steady
  uana = (beta \ (eye(ndof)-expm(-beta*nstep*Dt)))*(CP\S);
else
  uana = beta\(CP\S);
endif

Uana = real(phase*uana.');
                                #Uana=reshape(Uana,17,17)';
asave("uanalyt.tmp",Uana);
