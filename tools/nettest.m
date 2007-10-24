nodes = [1,
         2,
         3,
         4,
         5,
         6,
         8,
         9,
         10,
         11,
         12,
         13,
         14,
         15,
         16,
         17,
         18,
         19,
         20,
         21,
         22,
         23,
         24,
         25,
         26,
         27,
         28,
         29,
         30,
         31,
         32,
         33,
         34,
         35,
         36,
         37,
         38,
         39,
         40,
         41,
         42,
         43,
         44,
         45,
         46,
         47,
         48,
         49,
         50,
         51,
         52,
         65,
         66,
         67,
         69,
         70,
         71,
         72];

a=aloadg("nettest.log");

nprocs = max(a(:,2))+1;
bw = zeros(nprocs);

for row=a'
  j = row(1)+1;
  k = row(2)+1;
  bw(j,k) = row(3);
  bw(k,j) = row(3);
endfor

wmax = max(max(bw));
for k=1:nprocs
  bw(k,k)=wmax;
endfor

wset = [1:15,25:52]';

bwmin = 350;
format +
bw(wset,wset)>bwmin;
format 
all(bw(wset,wset)>bwmin) || warning("not correct wset");

for k=1:length(wset);
  printf("node%d\n",nodes(wset(k)));
endfor
