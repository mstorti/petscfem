L=20;
N=80;
x=L*((0:N)'/N*2-1);
X=x*ones(size(x))';
Y=X';
R=sqrt(X.^2+Y.^2);
n=10;
M = moviein(n);
for j=(0:n-1)
  fase=2*pi*j/n;
  surfl(x,x,sin(2*(R-fase))./(1+.3*R.^2))
  shading interp
  lighting gouraud
  lighting phong
  colormap(bone)
  axis([-L L -L L -1 40])
  M(:,j+1) = getframe;
end
movie(M,100,10)
       
