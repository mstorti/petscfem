#usage: 
function x = mapbouvt (nodi,nodj,xi,XNOD)

  if nodi>nodj
    x = mapbouvt(nodj,nodi,1-xi,XNOD);
    return;
  endif

  x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
endfunction
