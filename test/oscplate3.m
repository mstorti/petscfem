##__INSERT_LICENSE__
## $Id: oscplate3.m,v 1.3 2003/01/08 15:49:03 mstorti Exp $
source("~/.octaverc");

u16=aload("oscplate_16.sal");
u32=aload("oscplate_32.sal");
u128=aload("oscplate_128.sal");

printf("norm(u16-u128): %f\n",norm(u16-u128));
printf("norm(u32-u128): %f\n",norm(u32-u128));
printf("norm(u32-u128)/norm(u16-u128): %f\n",norm(u32-u128)/norm(u16-u128));

printf("norm(u32-u128)/norm(u16-u128) < 0.26 OK ? > %d \n",
       norm(u32-u128)/norm(u16-u128) < 0.26);
