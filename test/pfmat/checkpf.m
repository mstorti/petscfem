a=aload("data.tmp");
chkkey = a(1,1);
a(1,:)=[];

v = max(a);
v(1)==v(2) || error("matrix not square");

N=v(1)+1;
A=zeros(N);
for v=a';
  A(v(1)+1,v(2)+1)= A(v(1)+1,v(2)+1)+ v(3);
endfor

b=ones(N,1);
x = A\b;

x=[chkkey;
   x];
asave("xsol.tmp",x);
