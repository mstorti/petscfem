fid = fopen("nonlr.peri","w");
for k=1:2:35
  for j=1:5
    fprintf(fid,"%d %d %f   %d %d %f\n",k+1,j,-1.,k,j,1.);
    endfor
endfor
fclose(fid);
