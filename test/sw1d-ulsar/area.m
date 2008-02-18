B1=2;
B2=4;
Z1=1;
h=(0.:0.01:3)';
wl_width=[];
perimeter=[];
area=[];
for j=1:length(h)
  if (h(j)<=Z1)
    wl_width(j)  = B1;
    area(j)      = h(j)*B1;
    perimeter(j) = 2*h(j)+B1;
  else
    wl_width(j)  = B2;
    area(j)      = B1*Z1+B2*(h(j)-Z1);
    perimeter(j) = 2*h(j)+B2;
  endif
endfor
