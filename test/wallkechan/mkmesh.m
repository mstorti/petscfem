#__INSERT_LICENSE__

source("data.m.tmp");

hav=L/N;
w=zhomo([0 L 0 Ly],N+1,Ny+1,[1 0 hratio 1 0 hratioy]);
## last Ny_eq rows will be equal in size
wl=w(Ny+1+(-Ny_eq:0)',:);
xx=real(wl);
yy=imag(wl);
yy=yy(:,1);
eta=(0:Ny_eq)'/Ny_eq;
yy=yy(1)+eta.*(1+0.*(1-eta)).*(yy(Ny_eq+1)-yy(1));
yy=yy(:,ones(1,size(xx,2)));
wl=xx+i*yy;
w(Ny+1+(-Ny_eq:0)',:)=wl;
xx=xx(1,:)';

[xnod,icone]=pfcm2fem(w);

xnod = [xnod;
        xnod(1:Ny+1,:)];
ywall = xnod(1:Ny+1,2);
#ywall = (ywall(1:Ny)+ywall(2:Ny+1))/2;
if uwall_match_len>0 && piston
  xi = ywall/uwall_match_len;
  uwall = (3/4*((2*xi-1)-(2*xi-1).^3/3)+.5).*(xi<1) + (xi>=1);
else
  uwall = ones(size(ywall));
endif
epsuw = .1;
uwall=(1-epsuw)*(1-uwall);

icone=icone(:,[1 4 3 2]);
if peri
  nonlr=[(1:Ny)',(N+1)*(Ny+1)+(1:Ny)',0*uwall,uwall];
  lagr = (N+1)*(Ny+1)+(1:Ny);
else
  nonlr=[(1:Ny+1)',(N+1)*(Ny+1)+(1:Ny+1)',0*uwall,uwall];
  lagr = (N+1)*(Ny+1)+(1:Ny+1);
endif
## Don't include the first `nonlr' element, since all variables are
## fixed there 
#nonlr(rows(nonlr),:)=[]; 
nonlr(1,:)=[];
lagr(1)=[];
wallke = [(1:Ny)' (2:Ny+1)' 0*ones(Ny,1) xcent(uwall)];

fid = fopen("wallke.fixa_lag.tmp","w");
for k = lagr
  for j=3:5
    fprintf(fid,"%d %d %f\n",k,j,0.);
  endfor
endfor
fclose(fid);

epsp=0.;
uin=[0. Uav 0. 1e-3 1e-5];  
fid = fopen("wallke.fixa_in.tmp","w");
for jj=1:N+1
  k=1+(jj-1)*(Ny+1);
  for j=[1 4 5]
    fprintf(fid,"%d %d %f\n",k,j,uin(j));
  endfor
  if piston
    fprintf(fid,"%d %d %f\n",k,2,Uav);
  else
    fprintf(fid,"%d %d %f\n",k,2,Uav*(epsp+xx(jj))*(2*L+epsp-xx(jj))/L^2*1.5);
  endif
endfor
fclose(fid);

inlet = (1:(Ny+1):(N+1)*(Ny+1))';
outlet = (Ny+1:(Ny+1):(N+1)*(Ny+1))';

fid = fopen("wallke.peri.tmp","w");
fiduout = fopen("wallke.u_out_0.tmp","w");
fidpout = fopen("wallke.p_out_0.tmp","w");
fidp2out = fopen("wallke.wp_out_0.tmp","w");
fidfop = fopen("wallke.fixa_out_peri.tmp","w");
r=.5;
for k=1:N+2
  in = (k-1)*(Ny+1)+1;
  out = k*(Ny+1);
  fprintf(fiduout,"%d %d %f\n",out,1,0.);
  fprintf(fidpout,"%d %d %f\n",out,3,0.); 
  fprintf(fidp2out,"%f %d %d    %f %d %d\n",
          1,out,3,
          r,out-1,3,
          r^2,out-2,3);
  for j=[1 2 4 5]
    fprintf(fid,"%f %d %d    %f %d %d \n",-1.,out,j,+1.,in,j);
##    fprintf(fidfop,"%f %d %d    %f %d %d \n",-1.,out,j,+1.,out-1,j);
    fprintf(fidfop,"%f %d %d    %f %d %d    %f %d %d \n",-1.,out,j,+2.,out-1,j,-1.,out-2,j);
  endfor
endfor
fclose(fid);
fclose(fiduout);
fclose(fidpout);
fclose(fidp2out);
fclose(fidfop);

fid = fopen("wallke.wallke.tmp","w");
for k=1:rows(wallke)
  fprintf(fid,"%d %d %f %f\n",wallke(k,:));
endfor
fclose(fid);

fid = fopen("wallke.nonlr.tmp","w");
for k=1:rows(nonlr)
  fprintf(fid,"%d %d %f %f \n",nonlr(k,:));
endfor
fclose(fid);

uini=[0. Uav 0. uin(4:5)];
uini=[uini(ones((N+1)*(Ny+1),1),:);
      zeros(Ny+1,5)];

boux0=[(1:Ny+1)' ones(Ny+1,1) zeros(Ny+1,1)];
boux1=[((Ny+1)*N+(1:Ny+1))' ones(Ny+1,1) zeros(Ny+1,1)];
bovx0=[(1:Ny+1)' 2*ones(Ny+1,1) 0.*ones(Ny+1,1)];

asave("wallke.nod.tmp",xnod);
asave("wallke.con.tmp",icone);
asave("wallke.ini.tmp",uini);
#asave("wallke.nonlr.tmp",nonlr);
#asave("wallke.wallke.tmp",wallke);
asave("wallke.u0.tmp",boux0);
asave("wallke.u_0_right_wall.tmp",boux1);
asave("wallke.v_0_left_wall.tmp",bovx0);
