global parab_param
parab_param.a = 1.3;
parab_param.expo = 4;
parab_param.coe2 = 0.2;

x1 = [parabf(2) 2]';
x2 = [parabf(-2) -2]';

xi = (0:40)'/40;
x = zeros(rows(xi),2);

for j=1:rows(xi)
  xpro = project(xi(j),x1,x2,"parab");
  x(j,:) = xpro';
endfor
