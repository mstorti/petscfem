%clear all;
nod=load('canal1driem.nod.tmp');
nnod=size(nod,1)-2;
x=nod(1:nnod-2,1);
y=nod(1:nnod-2,2);
h_i=nod(1:nnod-2,3);
con=load('canal1driem.con.tmp');
%con=[con ones(size(con,1),1)];
dt=0.5;
field=2;Lx=20;

markersize_s=6;
linewidth_s=2.5;
fontsize_s=15;
fs=[1 1 1 10 10 10 10];

figure(1);
for i=0:1000
  i
%  fil=strcat({'STEPS-doble-rect-abso/canal1d.'},int2str(i),{'.tmp'});
%  fil2=strcat({'STEPS-doble-rect-noabso/canal1d.'},int2str(i),{'.tmp'});
  fil=strcat({'STEPS/canal1driem.'},int2str(i),{'.tmp'});
  u=load(char(fil));
  [ah]=axes;
%  [ah]=subplot(2,1,1);
  plot(x(1:nnod-2),u(1:nnod-2,field),'r','MarkerSize',markersize_s,'Linewidth',linewidth_s);
  grid;
  set(ah,'Fontsize', fontsize_s);
  set(ah,'XMinorGrid','off','YMinorGrid','off');
  xlabel('longitudinal coord x[m]');
  ylabel('free surf. height [m]');
  title('Riemann Absorbent B.C. for rectang channel shape');
  anot=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
  oh=text(1,1.035,anot);
  set(oh(1),'FontSize',fontsize_s);
  if (field==2)
    axis([0 Lx 0.98 1.45]);
  else
    axis([0 Lx -0.5 0.5]);
  end

%   [ah2]=subplot(2,1,2);
%   plot(x(1:nnod),unoabs(1:nnod,field),'b','MarkerSize',markersize_s,'Linewidth',linewidth_s);
%   grid;
%   set(ah2,'Fontsize', fontsize_s);
%   set(ah2,'XMinorGrid','off','YMinorGrid','off');
%   xlabel('longitudinal coord x[m]');
%   ylabel('free surf. height [m]');
%   title('Classical B.C. for circular channel shape');
%   anot2=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
%   oh2=text(1,1.035,anot2);
%   set(oh2(1),'FontSize',fontsize_s);
%   if (field==2)
%     axis([0 Lx 0.98 1.045]);
%   else
%     axis([0 Lx -0.5 0.5]);
%   end

%   fig=char(strcat({'YUV/.'},int2str(i),{'.tiff'}));
%   print('-dtiff',fig,'-zbuffer');
%   fig2=char(strcat({'YUV/doble-rect.'},int2str(i),{'.jpg'}));
%   print('-djpeg100',fig2,'-zbuffer');

pause;
hold off;
end
