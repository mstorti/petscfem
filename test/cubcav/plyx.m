%usage: 
function plyx (k,u)

  global N x3
  yy=[];
  uu=[];
  for ii=1:N+1
    indx = indxy(ii,k);
    yy = [yy x3(indx,2)];
    uu = [uu u(indx)];
  endfor
  plot(uu,yy);

endfunction
