#usage: 
function x = mapbou (nodi,nodj,xi,XNOD)

  global Rint Rext L Rmean

  if nodi>nodj
    x = mapbou (nodj,nodi,1-xi,XNOD);
    return;
  endif

  if nodi==1 && nodj==9
    theta = (1-xi)*pi/4;
    x = Rint*[cos(theta) sin(theta)];
  elseif nodi==1 && nodj==6
    theta = (1+2*xi)*pi/4;
    x = Rint*[cos(theta) sin(theta)];
  elseif nodi==2 && nodj==7
    theta = (1+2*xi)*pi/4;
    x = Rmean*[cos(theta) sin(theta)];
  elseif nodi==3 && nodj==8
    if xi<.5
      XX = [0 Rext];
      x = XNOD(3,:)+xi*(XX-XNOD(3,:));
    else
      theta = (1+2*xi)*pi/4;
      x = Rext*[cos(theta) sin(theta)];
    endif
  else
    x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
  endif
endfunction
