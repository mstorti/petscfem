#__ENDS_READING_VARS__

hratio=10;			# ratio of hmax/hmin
Ly=2;				# Lengh of mesh in y
slope=1e-4;			# bottom slope
h=0.1;				# water height
N=20;				# number of elements in y direction

source("../shallowt.m");
q=sqrt(Chezy^2*h*slope);
Fr=q/sqrt(gravity*h);

Pk=gravity*q^3/Chezy^2;
epsilon=Pk/h;

Pe=C_2*sqrt(C_mu)*gravity^(5/4)*q^4/(h*sqrt(D)*Chezy^2.5);
k=C_2*epsilon^2*h/Pe;

asave("tc.ini",kron([q*h 0 h h*k h*epsilon],ones(2*(N+1),1)));
tmpf="tc.data.tmp";
fid = fopen(tmpf,"w");
if !fid
  error(["Couldn't open " tmpf]);
endif

fprintf(fid,"$q=%f;\n",q);

fclose(fid);

##----<*>----<*>----<*>----<*>----<*>----<*>----<*>----<*>---- 
## Makes mesh
hav=Ly/N;
w=zhomo([0 Ly 0 hav],N+1,2,[1 0 hratio 1 0 1]);
[xnod,icone]=pfcm2fem(w);
#icone=icone(:,[1 4 3 2]);
xnod=xnod(:,[2 1]);
xnod=[xnod zeros(rows(xnod),1)];
fid = fopen("turbchanw.peri.tmp","w");
for k=1:N+1
  for kd=1:5
    fprintf(fid,"%f %d %d %f %d %d\n",-1,2*k,kd,+1,2*k-1,kd);
  endfor
endfor
fclose(fid);

asave("turbchanw.nod.tmp",xnod);
asave("turbchanw.con.tmp",icone);
