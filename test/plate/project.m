function xpro = project (xi,x1,x2,fun);

  xini = x1+xi*(x2-x1);
  nor = x2-x1;
  nor = [-nor(2); nor(1)];
  
  nor = nor/l2(nor);
  xpro = projectd(xini,nor,0.1,fun);

endfunction
