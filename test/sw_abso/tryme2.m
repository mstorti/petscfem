global g
g = 9.8;

u1=2;
h1=3;
wplus1 = u1+2*sqrt(g*h1);

u2=2.1;
h2 = ((wplus1-u2)/2)^2/g;
wplus2 = u2+2*sqrt(g*h2);

F1 = flux(u1,h1);
F2 = flux(u2,h2);
