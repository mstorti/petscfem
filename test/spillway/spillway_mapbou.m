#usage: 
function x = spillway_mapbou (nodi,nodj,xi,XNOD)

  global spillway_data

  if nodi>nodj
    x = parabolic_mapbou (nodj,nodi,1-xi,XNOD);
    return;
  endif

  if nodi==1 && nodj==2
    x = project(xi,XNOD(1,:)',XNOD(2,:)',"spillway_eq")';
  elseif nodi==4 && nodj==5
    x = project(xi,XNOD(4,:)',XNOD(5,:)',"fs_eq")';
  elseif nodi==5 && nodj==6
    x = project(xi,XNOD(5,:)',XNOD(6,:)',"fs_eq")';
  else
    x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
  endif
endfunction
