##__INSERT_LICENSE__
## $Id: check_fs.m,v 1.6 2005/09/25 19:44:51 mstorti Exp $
source("data.m.tmp");

load -force sqcav.ny.tmp
u = aload(save_file);
uref = aload(ref_file);
rem(N,2)==0 || error("N should be even!!");
ny = N/2*(N+1)+(1:N+1)';
u = u(ny,:);
erro = merr(u-uref);
tol=1e-10;
## plot([uref(:,1),u(:,1)],yh);

printf(["Square cavity at Re=1000. " \
	"Error < tol OK ? %d (error = %g, tol = %g)\n"], \
	erro<tol,erro,tol);
