X=aload("cubcav.nod.tmp");

x= 2*X(:,1)-1;
y= 2*X(:,2)-1;
z= 2*X(:,3)-1;

bubble = (1-x.^2).*(1-y.^2).*(1-z.^2);

v = sin(pi*x).*bubble;
w = sin(pi*y).*bubble;
u = sin(pi*x).*bubble;

asave("cubcav.state",[u v w 0*u]);
