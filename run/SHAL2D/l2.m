function r=l2(d)

  %if nargin==1
    r=zeros(size(d,1),1);
    r=sqrt(sum(d.^2));