nr=10;
nphi=30;
Rint=1;
Rext=5;

w=zhomo([log(Rint) log(Rext) -pi pi],nr+1,nphi+1,[1 0 5 1 0 1]);
w=-exp(w);
[xnod,icone]=pfcm2fem(w);

nnod=rows(xnod);
wake = [(1:nr+1:
