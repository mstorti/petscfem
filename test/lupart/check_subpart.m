##__INSERT_LICENSE__
## $Id: check_subpart.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
u_lu=aload("save.state.lu.tmp");
u_subpart=aload("save.state.sub_part.tmp");

tol = 1e-10;

err = merr(u_subpart-u_lu);
printf("IISD/Subpartitioning  OK ? > %d, [error: %g]\n",err<tol,err);
