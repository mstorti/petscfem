%usage: sA= rowlefs (s,A)
function sA= rowlefs (s,A)

  [N,m]=size(s);
  [NN,mp]=size(A);
  
  if N ~= NN
    error('row size of s and A should be the same')
  end

  p=round(mp/m);
  if mp ~= m*p
    error('rows(A) should be a multiple of rows(s)')
  end
  
  sA=zeros(size(A));
  for k=1:p
    sA(:,k:m:mp)=leftscal(s(:,k),A(:,k:m:mp));
  end
