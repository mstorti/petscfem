x = aload("step.nod.tmp");
dx = aload("step.state.tmp");
icone = aload("step.con.tmp");
s=1;
gplfem(x+s*dx,icone,"malla.gpl");
