function [uu,hh]= shall (u0,h0,DH)

  g=1;
  m=u0*h0;
  E0 = g*h0+0.5*u0^2;

  uu=roots([.5 0 g*DH-E0 g*m]);
  indx = find(uu<0);
  uu(indx) = [];
  hh = m./uu;
  
endfunction
