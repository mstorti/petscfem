## usage: 
## function Uri = primi2ri (Uprimi,gasdata)
function Uri = primi2ri (Uprimi,gasdata)

  ndof=columns(Uprimi);

  gamma = gasdata.gamma;
  Uri = zeros(size(Uprimi));

  if 0
    a = sqrt(gamma*Uprimi(:,4)./Uprimi(:,1));
    Uri(:,1) = Uprimi(:,2)-2*a/(gamma-1);
    Uri(:,2) = Uprimi(:,2)+2*a/(gamma-1);
    Uri(:,3) = log(Uprimi(:,4))-gamma*log(Uprimi(:,1));
    Uri(:,4) = Uprimi(:,3);
  else
    uref = gasdata.uref;
    rhoref = gasdata.rhoref;
    cref = gasdata.cref;
    pref = gasdata.pref;
    dU = Uprimi;
    dU(:,1) = dU(:,1) - rhoref;
    dU(:,2) = dU(:,2) - uref;
    dU(:,4) = dU(:,4) - pref;
    Uri(:,1) = dU(:,2) - dU(:,4)/(rhoref*cref);
    Uri(:,2) = dU(:,2) + dU(:,4)/(rhoref*cref);
    Uri(:,3) = dU(:,1) - dU(:,4)/cref^2;
    Uri(:,4) = dU(:,3);
  endif

endfunction
