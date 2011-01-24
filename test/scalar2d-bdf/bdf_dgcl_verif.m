###key bdf_dgcl_gaussian_diag_verif.m
### $Id: $
source("data.m.tmp");

u = aload(sprintf("gascont.state.%s.tmp",adv_case));

tol = 1e-7;
erro = merr(u-1);
printf("test OK? %d (erro %f, tol %f, \"adv_case %s\")\n",
       erro<tol,erro,tol,adv_case);
