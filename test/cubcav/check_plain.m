#$Id: check_plain.m,v 1.2 2003/12/08 13:19:00 mstorti Exp $

u0 = aload("cubcav.state.plain_bupl0.tmp");
u1 = aload("cubcav.state.plain_bupl1.tmp");

erro = merr(u1-u0);
tol = 1e-10;
printf("Cubcav/Block uploading error OK? > %d (error %g, tol %g)\n",
       erro<tol,erro,tol);

if 0
  ## Creates reference state file
  indx=sort(ceil(rand(30,1)*rows(u)));
  for k=indx'; 
    printf("%5d %20.12e %20.12e %20.12e %20.12e\n",k,u(k,:));
  endfor
  return
endif

uref = aload("cubcav.sol");
indx = uref(:,1);
uref = uref(:,2:5);

erro = merr(u0(indx,:)-uref);
printf("Cubcav/error wrt.reference solution OK? > %d (error %g, tol %g)\n",
       erro<tol,erro,tol);
