function fixag(f,dof,val,file)

  n=length(f);
  ff=[f dof*ones(n,1) val*ones(n,1)];
  asave(file,ff);

endfunction 
