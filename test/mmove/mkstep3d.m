source("data.m.tmp");
w = zhomo([0 1 0 1],N+1,N+1);
[xnod,icone] = pfcm2fem(w);
[xnod,icone] = extrude(xnod,icone,N,1/N);

#icone=[icone(:,[1 3 2]);
#       icone(:,[1 3 4])];

asave("step3d.nod.tmp",xnod);
asave("step3d.con.tmp",icone);

tol=1e-5;
x0 = find(abs(xnod(:,1))<tol);
x1 = find(abs(xnod(:,1)-1)<tol);
y0 = find(abs(xnod(:,2))<tol);
y1 = find(abs(xnod(:,2)-1)<tol);
z0 = find(abs(xnod(:,3))<tol);
z1 = find(abs(xnod(:,3)-1)<tol);

fid = fopen("step3d.fixa.tmp","w");
fixed = create_set([z0;z1;x0;x1;y0;y1]);

rmin = 0.5; rmax = 0.9;
for k=fixed
  x = xnod(k,1);
  y = xnod(k,2);
  z = xnod(k,3);
  r = sqrt(x^2+y^2);
  if r<rmin
    dd = 0;
  elseif r>rmax
    dd = -disp;
  else
    dd = -disp*(r-rmin)/(rmax-rmin);
  endif
  fprintf(fid,"%d %d %f\n",k,1,0.);
  fprintf(fid,"%d %d %f\n",k,2,0.);
  fprintf(fid,"%d %d %f\n",k,3,dd*z);
endfor
fclose(fid);

system("../../tools/hexasplit.bin -i step3d.con.tmp -o step3d.con-tet.tmp");
piecewtanh(2,"piecewise.dat.tmp");
