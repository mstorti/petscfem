#usage: 
function us = smsmooth (u,n_smoth_steps,omega_smooth);

N = rows(u)-1;
M = columns(u)-1;

u2 = u;
for j=1:n_smoth_steps
  u(:,M+1) = u(:,1);
  u(N+1,:) = u(1,:);

  for i=1:M+1;
    ip=modulo(i,N)+1;
    im=modulo(i-2,N)+1;
    u2(i,:) = omega_smooth*u(i,:)+(1-omega_smooth)/2*(u(ip,:)+u(im,:));
  endfor
  u = u2;

  for i=1:N+1;
    ip=modulo(i,M)+1;
    im=modulo(i-2,M)+1;
    u2(:,i) = omega_smooth*u(:,i)+(1-omega_smooth)/2*(u(:,ip)+u(:,im));
  endfor
  u = u2;

endfor
us = u;

endfunction
