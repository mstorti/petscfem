%clear all;
nod=load('plano.nod.tmp');
nnod=size(nod,1)-2;
x=nod(1:nnod,1);
y=nod(1:nnod,2);
h_i=nod(1:nnod,3);
con=load('plano.con.tmp');
%con=[con ones(size(con,1),1)];
dt=0.0133333333333333;
field=1;
g=1;
figure(1);
for i=0:1000
  i;
  fil=strcat({'STEPS/plano.'},int2str(i),{'.tmp'});
  u=load(char(fil));
  if (field==7)
    h1=u(1:nnod,1)./sqrt(g.*u(1:nnod,3));
    h2=u(1:nnod,1+3)./sqrt(g.*u(1:nnod,3+3));
  end
  if (field==3)
    h1=u(1:nnod,field);
    h2=h1+u(1:nnod,field+3);
  end
  if (field==1 | field==2)
    h1=u(1:nnod,field);
    h2=u(1:nnod,field+3);
  end
  plot(x(1:nnod),h1,'r-');
  hold on;
  plot(x(1:nnod),h2,'b-');
  if(field==1)
    axis([0 1 0.35 0.55]);
    grid;
    anot=char(strcat({'real time: '},num2str(i*dt,4),{' secs. '}));
    text(0.4,-0.65,anot);
    legend('fluido 1','fluido 2',3);
    title('velocities u1 and u2 vs x');
    xlabel('x [m]'); ylabel('u1 y u2 [m/sec]');
    box
    fig=char(strcat({'YUV/plano.'},int2str(i),{'.tiff'}));
    print('-dtiff',fig);
  end
  if(field==3)
    axis([0.45 0.9 0.95 1.05]);
    grid;
  end
  %  axis([0 2 0.8 2.8]);
%  pause;
  pause(0.01);
  hold off;
end
if (0)
  for i=0:1
    i
    fil=strcat({'STEPS/plano.'},int2str(i),{'.tmp'});
    u=load(char(fil));
    pdesurf([x y]', con', u(1:nnod,field));
				%  axis([0 1 0 0.20 0.1 0.9]);%view(2);
    view([0.7133    0.7009   -0.0000   -0.7071;
	  -0.4315    0.4391    0.7880   -0.3978;
	  -0.5523    0.5620   -0.6157    8.9632;
          0         0         0    1.0000]);
				%  caxis([0.41 0.55]);
    caxis([0. 0.15]);
    colormap('default');
    colorbar;
    anot=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
    text(0,0,0.3,anot);
    title('Hidraulic head vs time');
    xlabel('x [m]'); ylabel('y [m]'); zlabel('height [m]');
    box
    fig=char(strcat({'YUV/plano.'},int2str(i),{'.tiff'}));
    print('-dtiff',fig);
  end
end
