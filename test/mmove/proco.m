##__INSERT_LICENSE__
## $Id: proco.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
## x = aload("step.nod.tmp");
x = aload("lastmesh.dat");
## dx = aload("save.state.tmp");
icone = aload("step.con.tmp");
gplfem(x,icone,"malla.gpl");
