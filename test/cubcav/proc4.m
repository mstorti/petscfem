#usage: 
function  proc4 (u,N,x)
  np = (N+1)^2;
  minu = min(u);
  maxu = max(u);
  axis([0 1 minu maxu]);
  for k=2:N; 
    plot(x,reshape(u((k-1)*np+(1:np),1),N+1,N+1)); 
    pause; 
  endfor
  axis;
endfunction
