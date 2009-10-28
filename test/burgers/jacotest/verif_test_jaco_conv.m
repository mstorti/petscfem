###key verif_test_jaco_conv.m
### $Id: $

q = readconv("norm_res","test_jaco_conv_output.tmp");
printf("Convergence history:\n");
q

r1 = q(2)/q(1);
r2 = q(3)/q(2);
r3 = q(4)/q(3);

printf("convergence rate 1st to 2nd iter %g\n",r1);
printf("convergence rate 1st to 2nd iter %g\n",r2);
printf("convergence rate 1st to 2nd iter %g\n",r3);

Q1 = r2/r1^2;
Q2 = r3/r2^2;

printf("quadratic convergence test Q1 %g, OK %d\n",Q1,Q1<1);
printf("quadratic convergence test Q2 %g, OK %d\n",Q2,Q2<1);

printf("test OK %d\n",Q1<1 && Q2<1);
