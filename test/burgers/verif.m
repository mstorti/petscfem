##__INSERT_LICENSE__
## $Id: verif.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
data
CASE=getenv("CASE");
if strcmp(CASE,"_d001") 
  tol=1e-3;
elseif strcmp(CASE,"_d01") 
  tol=1e-4;
endif

uw=aload(["save.state." CASE "_wf.tmp"]);
unw=aload(["save.state." CASE "_nwf.tmp"]);
erro=merr(uw-unw);
printf("max. error %f, <tol(%f) OK? %d\n",erro,tol,erro<tol);
