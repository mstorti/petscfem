xnod = aload("cubcav.nod.tmp");
icone = aload("cubcav.con-tetra.tmp");

u = sum((xnod.^2)')';
x = rand(10,3);

ux = pfinterp(xnod, icone, u, x);

uxx = sum((x.^2)')';
[ux, uxx, ux-uxx]
