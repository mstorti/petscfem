function uL= supercr (uR)

  global gravity

  m=uR(1)*uR(2);
  UL=roots([.5 0 -(uR(2)^2/2+uR(1))  m]);
  UL=UL(find(UL>0));
  if length(UL) ~= 2
    error('Mas de dos soluciones positivas!!')
  end

  uL=zeros(size(uR));
  if abs(UL(1)-uR(2))>abs(UL(2)-uR(2))
    uL(2)=UL(1);
  else
    uL(2)=UL(2);
  end
    
  uL(1)=m/uL(2);
% FrL=uL(2)/sqrt(gravity*uL(1))
  uL(2)=m/uL(1);
