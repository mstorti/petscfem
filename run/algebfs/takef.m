function y = takef(x,indx,indx2);

  y=x(indx,:);
  x(indx,:)= 0*x(indx,:);

  if nargin>2
    x=y;
    y=x(:,indx2);
    x(:,indx2)= 0*x(:,indx2);
  end

  if max(abs(x))>0
    warning('Elementos no nulos');
  end

