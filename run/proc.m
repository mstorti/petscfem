U = aload("outvector.sal");
#U = aload("channel.ini");
U=U(1:102,:);

#xnod = aload("channel_stretch.nod");
xnod = aload("channel_fondo_var3.nod");
#xnod = aload("channel.nod");

nnod=rows(xnod);
x=xnod(1:2:102,1);
dh=xnod(1:2:102,3);

g = 1;
N=rows(U)/2-1;

U=U(1:2:2*N+2,:);

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

