%usage: 
function u = read_state (filename,nnod,ndof,indx)

  fid = fopen(filename,"r");
  if ! fid 
    error(["can't open " filename]);
  endif
  
  for k=1:indx
    [u,n] = fscanf (fid,"%lf",[ndof,nnod]);
  endfor
  fclose(fid);

  if n != nnod*ndof
    printf("read less than %d rows\n",nnod);
  endif
  u = u ';

endfunction
