## $Id: proc.m,v 1.1 2004/11/24 23:22:04 mstorti Exp $
a=[16 0 
   17 6 
   18 16 
   19 69 
   20 279 
   21 1001 
   22 4090 
   23 16282 ];
a(:,2) = a(:,2) +(a(:,2)==0);

x = log(a(:,1));
y = log(a(:,2));

A = [x, ones(size(x))];
param = (A'*A)\(A'*y);
aa = param(1);
bb = param(2);

yana = aa*x+bb;
## yana = 2.^(a*x+b);

loglog(2.^a(:,1),[a(:,2),e.^yana]);
## plot(a,[y yana])
