source("data.m.tmp");

phi=aload("save.state.tmp");
phi = reshape(phi,M+1,N+1);

xx=aload("advec.nod.tmp");
y=xx(1:M+1,2);
x=xx(1:(M+1):rows(xx),1);

phi_out = phi(:,N+1);
phi_anal = 2*(y>.75 )-1;
## plot(y,[phi_anal phi_out])

over_shoot = (max(max(phi_out),1)-1)/2
under_shoot = (max(-min(phi_out),1)-1)/2

## We accept a 5% over- and under- shoots
tol = 0.05;
#  if noise
#    tol = 0.1;
#  endif
printf("Over shoot < tol OK ? %d (over_shoot: %f, tol %f)\n",
       over_shoot<tol, over_shoot, tol);
printf("Under shoot < tol OK ? %d (under_shoot: %f, tol %f)\n",
       under_shoot<tol, under_shoot, tol);

## Number of nodes with an error > 5%
nodes_le = sum(abs(phi_out-phi_anal)>0.05*2);
max_num_ratio =  0.10;
max_num = ceil(max_num_ratio*(M+1));
printf("Number of nodes with error > tol OK ? %d\n",nodes_le<=max_num);
printf("node_le %d, max_num %d, max_num_ratio %f\n",
       nodes_le,max_num,max_num_ratio);

