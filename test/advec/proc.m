source("data.m.tmp");

x=aload("advec.nod.tmp");
phi_ini=aload("advec.ini.tmp");
[xc0,sx0,sy0]=bellpar(x,phi_ini);
phi=aload("save.state.tmp");
[xc,sx,sy]=bellpar(x,phi);

xrate = (sx/sx0-1)/(xc(1)-xc0(1))
yrate = (sy/sy0-1)/(xc(1)-xc0(1))
