##__INSERT_LICENSE__
## $Id: mkstep.m,v 1.9 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");
w = zhomo([0 1 0 1],N+1,N+1);
[xnod,icone] = pfcm2fem(w);

icone = icone(:,[1 4 3 2]);
if !strcmp(shape,"quad")
  icone=[icone(:,[2 3 1]);
	 icone(:,[4 1 3])];
endif

asave("step.nod.tmp",xnod);
asave("step.con.tmp",icone);

fid = fopen("step.fixa.tmp","w");
for k=1:N+1
  y = xnod(k,2);
  if y<.3
    dd = 0.;
  elseif y>.7
    dd = disp;
  else
    dd = disp*(y-.3)/.4;
  endif
  fprintf(fid,"%d %d %f\n",k,1,dd);
  fprintf(fid,"%d %d %f\n",k,2,0.);
  fprintf(fid,"%d %d %f\n",(N+1)*N+k,1,0.);
  fprintf(fid,"%d %d %f\n",(N+1)*N+k,2,0.);

  if k!=1 && k!=N+1
    node = (N+1)*(k-1)+1;
#    dd = disp * (1-xnod(node,1));
    dd = 0;
    fprintf(fid,"%d %d %f\n",node,1,dd);
    fprintf(fid,"%d %d %f\n",node,2,0.);
    node = (N+1)*k;

    dd = disp * (1-xnod(node,1));
    fprintf(fid,"%d %d %f\n",node,1,dd);
    fprintf(fid,"%d %d %f\n",node,2,0.);
  endif

endfor
fclose(fid);

piecewtanh(slope,"piecewise.dat.tmp");
