## $Id: cubcav_verif.m,v 1.1 2003/11/28 03:13:13 mstorti Exp $

u0 = aload("cubcav.state.plain_bupl0.tmp");
u1 = aload("cubcav.state.plain_bupl1.tmp");

merr(u1-u0)
