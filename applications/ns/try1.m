c1 = -5*log(5)+5;
c2 = 5*log(30)+c1-2.5*log(30);

f1=y;
f2=5*log(y)+c1;
f3=2.5*log(y)+c2;
f4=f1-ctff(f1-f3,3);
f5=wallf(y);
