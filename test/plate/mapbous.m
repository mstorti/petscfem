#usage: 
function x = mapbous(nodi,nodj,xi,XNOD)

  global Rint Rext L Rmean Rext2 

  if nodi>nodj
    x = mapbous (nodj,nodi,1-xi,XNOD);
    return;
  endif

  reflex = [1 -1];

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
      x = XNOD(3,:)+2*xi*(XX-XNOD(3,:));
    else
      theta = (1+2*xi)*pi/4;
      x = Rext*[cos(theta) sin(theta)];
    endif
  elseif nodi==6 && nodj==12
    x = mapbou_circle(6,12,xi,XNOD);
  elseif nodi==7 && nodj==13
    x = mapbou_circle(7,13,xi,XNOD);
  elseif nodi==8 && nodj==14
    x = mapbou_circle(8,14,xi,XNOD);
  else
    x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
  endif
endfunction
