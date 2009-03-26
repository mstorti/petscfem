clear all; close all;
nod=load('canal1d.nod.tmp');
nnod=size(nod,1)-2;
x=nod(1:nnod,1);
y=nod(1:nnod,2);
h_i=nod(1:nnod,3);
con=load('canal1d.con.tmp');
%con=[con ones(size(con,1),1)];
dt=0.5;
field=2;Lx=20;

markersize_s=6;
linewidth_s=2.5;
fontsize_s=15;
fs=[1 1 1 10 10 10 10];

for i=0:120
  i
  fil=strcat({'STEPS-doblerect-ulsar/canal1d.'},int2str(i),{'.tmp'});
  fil2=strcat({'STEPS-doblerect-noabso/canal1d.'},int2str(i),{'.tmp'});
%  fil=strcat({'STEPS-circ2-ulsar-new/canal1d.'},int2str(i),{'.tmp'});
%  fil2=strcat({'STEPS-circ2-noabso-new/canal1d.'},int2str(i),{'.tmp'});
%  fil=strcat({'STEPS-doblerect-ulsar/canal1d.'},int2str(i),{'.tmp'});
%  fil2=strcat({'STEPS-rect-riemann/canal1driem.'},int2str(i),{'.tmp'});
  uabs=load(char(fil));
  unoabs=load(char(fil2));
  figure(1);
% [ah]=axes;
  [ah]=subplot(2,1,1);
  plot(x(1:nnod),uabs(1:nnod,field),'r','MarkerSize',markersize_s,'Linewidth',linewidth_s);
  grid;
  set(ah,'Fontsize', fontsize_s);
  set(ah,'XMinorGrid','off','YMinorGrid','off');
  xlabel('longitudinal coord');
  ylabel('free surf. height');
  title('ULSAR B.C. for generic channel section');
  anot=char(strcat({'real time: '},num2str((i+1)*dt),{' secs. '}));
  oh=text(1,1.035,anot);
  set(oh(1),'FontSize',fontsize_s);
  if (field==2)
    axis([0 Lx 0.9 1.45]);
  else
    axis([0 Lx -0.5 0.5]);
  end

  [ah2]=subplot(2,1,2);
  plot(x(1:nnod),unoabs(1:nnod,field),'b','MarkerSize',markersize_s,'Linewidth',linewidth_s);
  grid;
  set(ah2,'Fontsize', fontsize_s);
  set(ah2,'XMinorGrid','off','YMinorGrid','off');
  xlabel('longitudinal coord');
  ylabel('free surf. height');
  title('Classical B.C. for generic channel section');
  anot2=char(strcat({'real time: '},num2str((i+1)*dt),{' secs. '}));
  oh2=text(1,1.035,anot2);
  set(oh2(1),'FontSize',fontsize_s);
  if (field==2)
    axis([0 Lx 0.9 1.45]);
  else
    axis([0 Lx -0.5 0.5]);
  end

%  fig=char(strcat({'YUV/circular.'},int2str(i),{'.tiff'}));
%  print('-dtiff',fig,'-zbuffer');
  fig2=char(strcat({'YUV/generic.'},int2str(i+1),{'.jpg'}));
  print('-djpeg100',fig2,'-zbuffer');

%  pause;
end
