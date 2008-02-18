clear all;
x=load('plano.nod.tmp');
u=load('STEPS/plano.10.tmp');
field=3;
plot3(x(:,1),x(:,2),u(:,field),'r.');

