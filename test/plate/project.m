## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: project.m,v 1.2 2003/03/14 14:32:41 mstorti Exp $

## usage: xpro = project (xi,x1,x2,fun);
## Finds a point on the curve defined by `fun(x)=0' in the interval
## delimited by `x1' and `x2', at a distance from both controlled by
## `xi'. More precisely, if `x0' is the projection of the point on the
## curve onto the secant from `x1' to `x2', then `xi=|x0-x1|/|x2-x1|'. 
## In particular, for `xi=0' we obtain `x1' and for `xi=1' we obtain `x2'.

## Author: Mario Storti
## Keywords: KEYWORDS
function xpro = project (xi,x1,x2,fun);

  xini = x1+xi*(x2-x1);
  nor = x2-x1;
  nor = [-nor(2); nor(1)];
  
  nor = nor/l2(nor);
  xpro = projectd(xini,nor,0.1,fun);

endfunction
