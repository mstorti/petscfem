x=(0:.1:10);
n=10;
M = moviein(n);
for j=(0:n-1)
  fase=2*pi*j/n;
  plot(x,sin(x+fase))
  M(:,j+1) = getframe;
end
movie(M)
       
