x0 = aload("step.nod.tmp");
dx = aload("step.state.tmp");
icone = aload("step.con.tmp");
x = x0+dx;
gplfem(x,icone,"malla.gpl");
