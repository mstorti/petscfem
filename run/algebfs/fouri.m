AA=aload("a.dat");
AA=getblock(AA);
ndof=3;

nx=sqrt(n);
nx2=nx/2;
nodo=nx*(nx2-1)+nx2;

ID=zeros(nx,nx);
ID(nx2,nx2)=1;
tID=fft2(ID);

At=zeros(n,ndof^2);
Ant=zeros(n,ndof^2);
cou=0;
for k=1:ndof
  for l=1:ndof
    cou=cou+1;
    akl=AA((nodo-1)*ndof+k,l:ndof:ndof*n);
    Ant(:,cou) = vec(akl);
    akl=reshape(akl,nx,nx);
    akl=fft2(akl)./tID;
    At(:,cou) = vec(akl);
  end
end

H=zeros(n,1);
Hprec=zeros(n,1);
gd=zeros(n,1);
giAd=zeros(n,1);

for k=1:n
  
  A=reshape(At(k,:),ndof,ndof).';
  AAA = A(1:2,1:2);
  B=A(1:2,3);
  C=A(3,1:2);
  DD=A(3,3);
  H(k) = -DD +C*inv(AAA)*B;
  Hprec(k) = (-DD +C*inv(AAA)*B);
  gd(k) = C*B;
  giAd(k) = C*inv(AAA)*B;
  keyboard 
endfor

H=reshape(H,nx,nx);
Hprec=reshape(Hprec,nx,nx);
gd=reshape(gd,nx,nx);
giAd=reshape(giAd,nx,nx);
D=reshape(At(:,9),nx,nx);
