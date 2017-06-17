###key proc2.m
### $Id: $

xnod = aload("./sqcav.nod.tmp");
icone = aload("./sqcav.con.tmp");
u = aload("./sqcav.weak_form_1.tmp");

xdmfsave("sqcav",xnod,icone,u,name);
