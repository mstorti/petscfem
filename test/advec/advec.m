source("data.m.tmp");

w=zhomo([0 1 0 1],N+1,M+1);
[xnod,icone]=pfcm2fem(w);

icone=[icone(:,[1 4 3]);
       icone(:,[3 2 1])];
	     
asave("advec.nod.tmp",xnod);
asave("advec.con.tmp",icone);

fid = fopen("advec.fixa.tmp","w");
for k=1:M+1
  y = xnod(k,2);
  phi = tanh((y-0.5)/delta);
  fprintf(fid,"%d %d %f\n",k,1,phi);
endfor
fclose(fid);
