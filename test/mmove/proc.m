x = aload("step.nod.tmp");
dx = aload("save.state.tmp");
icone = aload("step.con.tmp");
s=1;
gplfem(x+s*dx,icone,"malla.gpl");
