## usage: 
function diffujac(fid,len)

  ## $Id: diffujac.m,v 1.2 2007/02/20 15:08:26 mstorti Exp $
  fprintf(fid,"diffusive_jacobians_mol ");
  for k=1:len
    fprintf(fid," 0.0");
  endfor
  fprintf(fid,"\n");

endfunction
