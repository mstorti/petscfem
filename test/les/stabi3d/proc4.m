if 1
  dudy = uav(2,:)'/y(2);
  tauw = nu*dudy;
  ustar = sqrt(tauw);
  nsteps = columns(uav);
  U_bulk = zeros(nsteps,1);
  for k=1:nsteps
    U_bulk(k) = sum(xcent(uav(:,k)).*diff(y));
  endfor
  Cf = tauw./(0.5*U_bulk.^2);

else

  uav = sum(u(1:Ny+1,start:last)')'/nsteps;
  uav = sum(u(1:Ny+1,start:last)')'/(last-start+1);
  
  dudy = uav(2)/y(2);
  tauw = nu*dudy;
  ustar = sqrt(tauw);
  
  U_bulk = sum(xcent(uav).*diff(y));
  
  Cf = tauw/(0.5*U_bulk^2);

endif
