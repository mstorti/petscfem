#usage: 
function x = mapbou (nodi,nodj,xi,XNOD)

  global Rint Rext L Rmean Rext2 

  if nodi>nodj
    x = mapbou (nodj,nodi,1-xi,XNOD);
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
  elseif nodi==14 && nodj==19
    x = mapbou(3,8,xi,XNOD) .* reflex;
  elseif nodi==12 && nodj==17
    x = mapbou_circle(12,17,xi,XNOD);
  elseif nodi==13 && nodj==18
    x = mapbou_circle(13,18,xi,XNOD);
  elseif nodi==6 && nodj==17
    x = mapbou_circle(6,17,xi,XNOD);
  elseif nodi==7 && nodj==18
    x = mapbou_circle(7,18,xi,XNOD);
  elseif nodi==8 && nodj==19
    x = mapbou_circle(8,19,xi,XNOD);
  elseif nodi==21 && nodj==22
    if xi<0.5
      XX = [0 Rext2];
      x = XNOD(21,:)+2*xi*(XX-XNOD(21,:));
    else
      theta = (1+2*xi)*pi/4;
      x = Rext2*[cos(theta) sin(theta)];
    endif
  elseif nodi==22 && nodj==23
    x = mapbou_circle(22,23,xi,XNOD);
  elseif nodi==23 && nodj==24
    x = mapbou(21,22,1-xi,XNOD) .* reflex;
  else
    x = XNOD(nodi,:)+xi*(XNOD(nodj,:)-XNOD(nodi,:));
  endif
endfunction
