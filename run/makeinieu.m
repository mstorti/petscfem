xnod = aload("euler2d.nod");
x = xnod(:,1);
y = xnod(:,2);
r=sqrt(x.^2+y.^2);

N=861;
DP = 1;
sigma=0.5;

u = ones(N,1) * [1 1 0 180] ;
u(:,4) = u(:,4) + DP*exp(-(r/sigma).^2);
asave("euler2d.ini",u);
