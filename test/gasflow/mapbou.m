#usage: 
function x = mapbou (nodi,nodj,xi,XNOD)

  global Rin Rn nw L

  if nodi>nodj
    x = mapbou (nodj,nodi,1-xi,XNOD);
    return;
  endif

  if nodi==3 && nodj==4
    x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
    coef = acosh(10)/nw;
    x(2) = x(2) - (Rin-Rn)*(1/cosh(coef*x(1)) - 1/cosh(coef*L));
  else
    x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
  endif
endfunction
