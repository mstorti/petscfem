#usage: [amp,omega,phase]= trigreg (x,dt)
#   Find amp, phase so that (in a least squares sense)
#         x[j] = amp * sin (omega*(j-1)*dt+phase)
#
#    args:  x, the vector of x values
#           dt, time step (default=1)
#         
function [amp,omega,phase]= trigreg (x,dt)

  if nargin==1
    dt=1;
  endif

  N=rows(x);
  if N<3
    error("at least 3 values are required");
  endif

  [bid,j]=max(abs(x));
  d2x = (x(j+1)-2*x(j)+x(j-1))/dt^2;
  omega = sqrt(-d2x/x(j));

endfunction
