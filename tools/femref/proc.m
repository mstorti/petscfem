## $Id: proc.m,v 1.4 2004/11/25 00:11:47 mstorti Exp $

## N | collisions with hash | colls. with random
a=[16 0            2 	   
   17 6 	   6 	   
   18 16 	   12 	   
   19 69 	   63 	   
   20 279 	   238    
   21 1001 	   978    
   22 4090 	   4005   
   23 16282 	   16412  
   24 64843	   65652  
   25 261172	   260751 
   26 1037443      1037578];

a(:,2) = a(:,2) +(a(:,2)==0);

x = a(:,1)*log(2);
y = log(a(:,2));

A = [x, ones(size(x))];
param = (A'*A)\(A'*y);
aa = param(1);
bb = exp(param(2))*2^32;

yana = bb.*2^(-32)*(2.^a(:,1)).^aa;
## yana = 2.^(a*x+b);

loglog(2.^a(:,1),[a(:,2:3),yana]);
## plot(a,[y yana])

## produces aa=1.97, bb=1.52
## nbr of collisions is bb*m^aa/N
