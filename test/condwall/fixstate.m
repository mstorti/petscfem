## $Id: fixstate.m,v 1.1 2005/03/28 18:19:24 mstorti Exp $
source("data.m.tmp");

dx_step=str2num(getenv("dx_step"));
printf("dx_step %d\n",dx_step);

u = aload(sprintf("STEPS/condwall.state_%d.tmp",dx_step));

nnod = (Nx1+Nx2+2)*(Ny+1);
asave("condwall.dx-state.tmp",u(1:nnod,:));
