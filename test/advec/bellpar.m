#usage: 
function [xc,sx,sy]= bellpar (x,phi)

  xc = [sum(phi.*x(:,1)) sum(phi.*x(:,2))]/sum(phi);
  sx = sum(phi.*(x(:,1)-xc(1)).^2)/sum(phi);
  sy = sum(phi.*(x(:,2)-xc(2)).^2)/sum(phi);

endfunction
