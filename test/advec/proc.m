source("data.m.tmp");

x=aload("advec.nod.tmp");
phi_ini=aload("advec.ini.tmp");
[xc0,sx0,sy0]=bellpar(x,phi_ini);
printf("x0: %f %f, sx0 %f, sy0 %f\n",xc0(1),xc0(2),sx0,sy0);

phi=aload("save.state.tmp");
[xc,sx,sy]=bellpar(x,phi);
printf("x: %f %f, sx %f, sy %f\n",xc(1),xc(2),sx,sy);

xrate = (sx/sx0-1)/(xc(1)-xc0(1))
yrate = (sy/sy0-1)/(xc(1)-xc0(1))
