#usage: 
function uu = refin (u,ref)

  N = round(sqrt(rows(u)))-1;
  u = reshape(u,N+1,N+1);
  NN = ref*N;

  uu=zeros(NN+1);
  uu(1:ref:NN+1,1:ref:NN+1) = u;

  for k=2:NN
    kk = floor((k-1)/ref)*ref+1;
    alpha = (k-kk)/ref;
    uu(k,:) = (1-alpha)*uu(kk,:) + alpha * uu(kk+ref,:);
  endfor

  for k=2:NN
    kk = floor((k-1)/ref)*ref+1;
    alpha = (k-kk)/ref;
    uu(:,k) = (1-alpha)*uu(:,kk) + alpha * uu(:,kk+ref);
  endfor

  uu = smsmooth(uu,4,0.5);
  uu = vec(uu);

endfunction
