##__INSERT_LICENSE__
## $Id: sine_crank_nic.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
u16 = aload("save.state.16");
u32 = aload("save.state.32");
u128 = aload("save.state.128");

e16=merr(u16-u128);
e32=merr(u32-u128);

printf("||u_16-u_128|| = %f \n",e16);
printf("||u_32-u_128|| = %f \n",e32);
printf("||u_32-u_128|| / ||u_16-u_128|| = %f, < 0.25 OK? %d \n",e32/e16,e32/e16<0.25);
