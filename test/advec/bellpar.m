#usage: 
function [xc,sx,sy]= bellpar(x,icone,phi)

  sum_phi = pfint(x,icone,phi);
  xc = pfint(x,icone,leftscal(phi,x))/sum_phi;

  sx = pfint(x,icone,phi.*(x(:,1)-xc(1)).^2)/sum_phi;
  sy = pfint(x,icone,phi.*(x(:,2)-xc(2)).^2)/sum_phi;
endfunction
