## $Id: fixstate.m,v 1.2 2005/03/29 21:05:26 mstorti Exp $
source("data.m.tmp");

dx_step=str2num(getenv("dx_step"));
printf("dx_step %d\n",dx_step);

## sdir = "STEPS-2005-MAR-29";
sdir = "STEPS";
u = aload([sdir sprintf("/condwall.state_%d.tmp",dx_step)]);

nnod = (Nx1+Nx2+2)*(Ny+1);
asave("condwall.dx-state.tmp",u(1:nnod,:));
