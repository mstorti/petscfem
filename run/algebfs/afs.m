function rp = afs(p)

  global A B C D bu bp sU sp omega r0 homog_version

  p = p.*sp;

  but = bu - sparmul(B,p);
  u = zeros(size(bu));
  
  itmom=10;
  disp("Internal convergence: ");
  for k=1:itmom
    ru = but - sparmul(A,u);
    norm (ru)
    u = u + omega * (sU .* ru);
  endfor
  disp("Fin.");

  rp = sp.*(bp-sparmul(C,u)-sparmul(D,p));

  if homog_version
    rp=rp-r0;
  endif
  
endfunction 
