global AA A B C D b H bH

AA=aload("a.dat");
AA=getblock(AA);
N=rows(AA);
A2=AA(3:3:N,:);
A1=AA;

A1(3:3:N,:)=[];

D=A2(:,3:3:N);
A2(:,3:3:N)=[];
C=A2;
clear A2;

B=A1(:,3:3:N);
A1(:,3:3:N)=[];
A=A1;

clear A1
AA=[A B; 
    C D];

n=rows(D);
N=3*n;
b = rand(3*n,1);
b(1:2*n)=zeros(2*n,1);
II=2*n+(1:n)';
b(II)=b(II)-sum(b(II))/n;
bp=b(II);
b(N)=[];
uex = [AA(1:3*n-1,:);
       zeros(1,2*n) ones(1,n)] \ [b; 0];
uex = uex(1:N-1);

bu=b(1:2*n);

iA = inv(A);
H = C*iA*B-D;
bH = C*iA*bu-bp;
bHn=bH(n);
bH(n)=[];
pex = [H(1:n-1,:);ones(1,n)] \ [bH;0];

norm(pex(1:n-1)-uex(2*n+(1:n-1)))
