#usage: 
function x = mapbou (nodi,nodj,xi,XNOD)

  global Rint Rext L Rmean Rext2 

  x8 = XNOD(8,:)';
  x3 = XNOD(3,:)';
  x5 = XNOD(5,:)';
  x19 = XNOD(19,:)';
  x14 = XNOD(14,:)';
  x16 = XNOD(16,:)';

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
    x = project(xi,x3,x8,"parab")';
  elseif nodi==3 && nodj==5
    x = project(xi,x3,x5,"parab")';
  elseif nodi==9 && nodj==12
    x = mapbou_circle(9,12,xi,XNOD);
  elseif nodi==14 && nodj==19
    x = project(xi,x14,x19,"parab")';
  elseif nodi==14 && nodj==16
    x = project(xi,x14,x16,"parab")';
  elseif nodi==12 && nodj==17
    x = mapbou_circle(12,17,xi,XNOD);
  elseif nodi==13 && nodj==18
    x = mapbou_circle(13,18,xi,XNOD);
  elseif nodi==6 && nodj==17
    x = mapbou_circle(6,17,xi,XNOD);
  elseif nodi==7 && nodj==18
    x = mapbou_circle(7,18,xi,XNOD);
  elseif nodi==8 && nodj==19
    x = project(xi,x8,x19,"parab")';
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
