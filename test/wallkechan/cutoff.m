## usage:  w = cutoff (v,tol)
##    Sets components of `v' to `tol' if `v < tol',
##    else let unchanged. 

# $Id: cutoff.m,v 1.1 2001/07/01 15:14:26 mstorti Exp $
#__INSERT_LICENSE__

function v = cutoff (v,tol)

  v = tol + (v-tol).*((v-tol)>0);

endfunction
