n=round(T/dt/nsave);
if n>nsavemax
  error('Guarda que guardamos menos pasos de los que entran en un periodo!!')
end

% surfl contour plot
tipo='plot'; % 

clear M
M = moviein(n);
scale=10;
hinf=1;
%# v=scale*(minh+(*(maxh-minh)-hinf);

minh=min(min(hh));
maxh=max(max(hh));
dh=max([maxh-hinf hinf-minh]);

ncont=20; stretch=.99;
xi=2*(0:ncont)'/ncont-1;
v=hinf+dh*atanh(xi*stretch)/atanh(stretch);

for j=(0:n-1)
  x=(0:Nx)'*dx;
  y=(0:Ny)'*dy;
  jj=isave-j;
  if jj<=0
    jj=jj+nsavemax;
  end
  hhh=reshape(hh(:,jj),Nx+1,Ny+1)';
  if strcmp(tipo,'surfl')
    hhh=10*(hhh-1); % para escalear
    axis([0 Lx 0 Ly -1 1])
    axis manual
    surfl(x,y,hhh)
    shading interp
    lighting gouraud
    %lighting phong
    colormap(pink) % bone gray
    axis equal
  elseif strcmp(tipo,'plot')
    plot(x,hhh(1,:))
    axis([0 Lx minh maxh])
  elseif strcmp(tipo,'contour')
    contour(x,y,hhh,v)
    axis equal
    colormap(white)
    whitebg([0 0 0])
  else
    error('No existe ese tipo de ploteo!')
  end
  M(:,n-j) = getframe;
end
movie(M,100,10)
       
