##__INSERT_LICENSE__
## $Id: check_direct_superlu.m,v 1.3 2003/01/08 15:49:04 mstorti Exp $
if ! strcmp(getenv("SUPERLU"),"")
  u_lu=aload("save.state.lu.tmp");
  u_slu=aload("save.state.direct_superlu.tmp");
  
  tol = 1e-8;

  err = merr(u_slu-u_lu);
  printf("Direct/SuperLU  OK ? > %d, [error: %g]\n",err<tol,err);
else
  printf("Direct/SuperLU  OK ? > 1 (No SuperLU installed!)\n");
endif  
