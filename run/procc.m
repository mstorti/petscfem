U = aload("outvector.sal");
#U = aload("channel.ini");

indx = 1:21:861;
U=U(indx,:);

#xnod = aload("channel_stretch.nod");
xnod = aload("cholo.nod");
nnod=rows(xnod);
x=xnod(indx,1);
dh=xnod(indx,3);

g = 1;
N=rows(U)/2-1;

H=U(:,3);
U=U(:,1)./H;

title("h");
plot(x,H);
pause;

title("u");
plot(x,U);
pause

title("Fr");
Fr = U./sqrt(g*H);
plot(x,Fr);
pause

title("u*h");
plot(x,U.*H);
pause

title("g*(h+H)+0.5*u^2");
plot(x,g*(H+dh)+0.5*U.^2);

