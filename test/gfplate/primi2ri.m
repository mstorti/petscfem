## usage: 
## function Uri = primi2ri (Uprimi,gasdata)
function Uri = primi2ri (Uprimi,gasdata)

  ndof=columns(Uprimi);

  gamma = gasdata.gamma;
  a = sqrt(gamma*Uprimi(:,4)./Uprimi(:,1));
  Uri = zeros(size(Uprimi));
  Uri(:,1) = Uprimi(:,2)-2*a/(gamma-1);
  Uri(:,2) = Uprimi(:,2)+2*a/(gamma-1);
  Uri(:,3) = log(Uprimi(:,4))-gamma*log(Uprimi(:,1));
  Uri(:,4) = Uprimi(:,3);

endfunction
