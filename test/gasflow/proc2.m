q=aload("vtube.mass_flow_rate.tmp"); 
q=[q -.5*q(:,1)];
return

mtime=0; ff="";
for f=glob("nohup.out.*")'; 
  [info, erro, msg] = stat (f');
  if info.mtime>mtime
    mtime = info.mtime;
    ff = f';
  endif
endfor

ph = readconv("p0",ff);
pc = readconv("p1",ff);
hc = readconv("delta_u",ff);

n = min([length(ph) length(pc)]);
p=[ph(1:n) pc(1:n)];

title("Flow rates");
plot(q);
pause
title("Control pressures");
plot(p);
pause
title("Convergence history");
semilogy(hc);
title("");
