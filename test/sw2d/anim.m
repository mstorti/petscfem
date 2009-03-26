clear all;
nod=load('plano.nod.tmp');
nnod=size(nod,1);
x=nod(1:nnod,1);
y=nod(1:nnod,2);
h_i=nod(1:nnod,3);
%con=load('../mesh_slope_inicial/testsw2d.con');
con=load('plano.con.tmp');
con=[con ones(size(con,1),1)];
%dt=0.5;
minx=min(x);
maxx=max(x);
miny=min(y);
maxy=max(y);
field=3;
dt=0.15;
B_out = 0.25;
first=0;st=1;last=50;%1550;
for i=first:last
  i
  fil=strcat({'STEPS/plano.'},int2str(i),{'.tmp'});
%  fil=strcat({'STEPS-noslip-ulsaruref/testsw2d.'},int2str(i),{'.tmp'});
  u=load(char(fil));
  u=u(1:nnod,1:3);
  if (field==4)
    velmod=sqrt(u(1:nnod,1).^2+u(1:nnod,2).^2);
    u(:,field)=velmod;
  end
  if (field==3)
    u(:,field)=u(1:nnod,3)+h_i;
  end
  pdesurf([x y]', con', u(:,field));
%  pdesurf([x y]', con', velmod);
%  axis([minx maxx miny maxy]);% 0.4 0.6]);%view(2);
%  view([0.7133    0.7009   -0.0000   -0.7071;
%	-0.4315    0.4391    0.7880   -0.3978;
%	-0.5523    0.5620   -0.6157    8.9632;
%        0         0         0    1.0000]);
%  caxis([0.41 0.55]);
  caxis([0.9 1.5]);
  axis([0.0 1.0 0.0 1.0 0.50 1.5]);% 0.4 0.6]);%view(2);
%  axis equal;
  colormap('default');
  colorbar;
  anot=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
  text(0,0,0.3,anot);
  title('Hidraulic head vs time');
  xlabel('x [m]'); ylabel('y [m]');% zlabel('height [m]');
  box
  fig=char(strcat({'YUV/plano.'},int2str(i),{'.tiff'}));
  print('-dtiff',fig);
%  pause;
end
