##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.3 2003/01/08 15:49:04 mstorti Exp $
fid = fopen("nonlr.peri","w");
for k=1:2:35
  for j=1:5
    fprintf(fid,"%f %d %d   %f %d %d\n",-1.,k+1,j,1.,k,j);
    endfor
endfor
fclose(fid);
