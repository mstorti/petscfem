n=10;

Q=zeros(n,n);

Q(1:2,1:3) = ones(2,3);
Q(3:5,4:5) = ones(3,2);
Q(6:8,6:8) = eye(3);

[bid,II] = sort(rand(n,1));
[bid,JJ] = sort(rand(n,1));

Q=Q(II,JJ);
