data

u=aload("save.state.tmp");
if use_bcconv
  meu=merr(u-rval);
else 
  ua=aload("uanalyt.tmp");
  disp("max error over dof's")
  if fix_some_node
    for j=1:ndof
      u(:,j)=u(:,j)+ua(1,j);
    endfor
  endif
  meu=merr(abs(u-ua));
  disp(meu);
  us=reshape(u(:,1),ny+1,nx+1)';
  uas=reshape(ua(:,1),ny+1,nx+1)';
endif  

mu=max(max(abs(u)));
printf("max error over all dof's: %f\n",meu);
printf("rel. max error = max(err)/max(u) %g\n",meu/mu);
printf("rel. max error  <tol(%f) OK? %d\n",tol,meu/mu<tol);

