#usage: 
function x = mapbou_circle (nod1,nod2,xi,XNOD);

  r1 = l2(XNOD(nod1,:));
  r2 = l2(XNOD(nod2,:));
  abs(r1-r2)<1e-10 || error("Not on the same circunference");
  theta1 = atan2(XNOD(nod1,2),XNOD(nod1,1));
  theta2 = atan2(XNOD(nod2,2),XNOD(nod2,1));
  if abs(theta2-theta1)>pi
    theta2 = theta2-2*pi*sign(theta2-theta1);
  endif
  theta = theta1+xi*(theta2-theta1);
  x = r1*[cos(theta) sin(theta)];

endfunction
