data
u=aload("save.state.tmp");
ua=aload("uanalyt.tmp");
disp("max error over dof's")
if fix_some_node
  for j=1:ndof
    u(:,j)=u(:,j)+ua(1,j);
  endfor
endif
meu=max(abs(u-ua));
disp(meu);
printf("max error over all dof's: %f\n",max(meu));
mu=max(max(abs(ua)));
printf("rel. max error = max(err)/max(u) %g\n",max(meu)/mu);
printf("rel. max error  <tol(%f) OK? %d\n",tol,max(meu)/mu<tol);

us=reshape(u(:,1),ny+1,nx+1)';
uas=reshape(ua(:,1),ny+1,nx+1)';
