q=aload("vtube.mass_flow_rate.tmp"); 
q=[q -.5*q(:,1)];

ph = readconv("p0","nohup.out.8");
pc = readconv("p1","nohup.out.8");
hc = readconv("delta_u","nohup.out.8");

p=[ph pc];

title("Flow rates");
plot(q);
pause
title("Control pressures");
plot(p);
pause
title("Convergence history");
semilogy(hc);
title("");
