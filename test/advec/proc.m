source("data.m.tmp");

x=aload("advec.nod.tmp");
icone=aload("advec.con.tmp");
icone = icone(:,1:3);		# discard u columns in rotating cone example
phi_ini=aload("advec.ini.tmp");
## [xc0,sx0,sy0]=bellpar2(x,phi_ini);
[xc0,sx0,sy0]=bellpar(x,icone,phi_ini);
printf("x0: %f %f, sx0 %f, sy0 %f\n",xc0(1),xc0(2),sx0,sy0);

phi=aload("save.state.tmp");
## [xc,sx,sy]=bellpar2(x,phi);
[xc,sx,sy]=bellpar(x,icone,phi);
printf("x: %f %f, sx %f, sy %f\n",xc(1),xc(2),sx,sy);

