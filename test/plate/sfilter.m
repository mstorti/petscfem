omega=1;
dt=0.1;
xi=1;

a = [1+xi*omega+omega^2/2;
     -2;
     1-xi*omega+omega^2/2];

b = [1 0 1]';

x=ones(32,1);
x(1:10)=0;
y = filter(b,a,x);
