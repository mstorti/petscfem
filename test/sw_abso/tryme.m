g=9.8;
h=2;
u=0.5;

c = sqrt(g*h);

A=[2*u -u^2+g*h; 1 0];
eig(A)
u+c
u-c

aA = [2*(u^2+g*h),(g*h-u^2)*2*u; 2*u, 2*(g*h-u^2)]/(2*c);
eig(aA);


