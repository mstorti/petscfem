clear all; close all;
nod=load('canal1d.nod.tmp');
nnod=size(nod,1)-2;
x=nod(1:nnod,1);
y=nod(1:nnod,2);
h_i=nod(1:nnod,3);
con=load('canal1d.con.tmp');
%con=[con ones(size(con,1),1)];
dt=0.5;
field=1;Lx=20;

markersize_s=6;
linewidth_s=2.5;
fontsize_s=15;
fs=[1 1 1 10 10 10 10];

if (0)
  for i=0:1000
    i
    fil=strcat({'STEPS-rect-ulsar/canal1d.'},int2str(i),{'.tmp'});
    fil2=strcat({'STEPS-rect-riemann/canal1driem.'},int2str(i),{'.tmp'});
    uabs=load(char(fil));
    uriem=load(char(fil2));
    figure(1);
% [ah]=axes;
    [ah]=subplot(2,1,1);
    plot(x(1:nnod),uabs(1:nnod,field),'r','MarkerSize',markersize_s,'Linewidth',linewidth_s);
    grid;
    set(ah,'Fontsize', fontsize_s);
    %set(ah,'XMinorGrid','off','YMinorGrid','off');
    xlabel('longitudinal coord x[m]');
    ylabel('free surf. height [m]');
    title('ULSAR Absorbent B.C. for rectang. channel shape');
      anot=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
      oh=text(1,1.035,anot);
      set(oh(1),'FontSize',fontsize_s);
      if (field==2)
	axis([0 Lx 0.98 1.45]);
      else
	axis([0 Lx -0.5 0.5]);
      end
      
      [ah2]=subplot(2,1,2);
      plot(x(1:nnod),uriem(1:nnod,field),'b','MarkerSize',markersize_s,'Linewidth',linewidth_s);
      grid;
      set(ah2,'Fontsize', fontsize_s);
      set(ah2,'XMinorGrid','off','YMinorGrid','off');
      xlabel('longitudinal coord x[m]');
      ylabel('free surf. height [m]');
      title('Riemann Absorbent B.C. for rectang. channel shape');
	anot2=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
	oh2=text(1,1.035,anot2);
	set(oh2(1),'FontSize',fontsize_s);
	if (field==2)
	  axis([0 Lx 0.98 1.45]);
	else
	  axis([0 Lx -0.5 0.5]);
	end
	
	fig=char(strcat({'YUV/riem-comp.'},int2str(i),{'.tiff'}));
	print('-dtiff',fig,'-zbuffer');
	fig2=char(strcat({'YUV/riem-comp.'},int2str(i),{'.jpg'}));
	print('-djpeg100',fig2,'-zbuffer');
	
%  pause;
  end
end

ut=[];uriemt=[];
if (0)
  for i=0:299
    i
    fil=strcat({'STEPS-rect-nouref/canal1d.'},int2str(i),{'.tmp'});
    fil2=strcat({'STEPS-rect-riemann/canal1driem.'},int2str(i),{'.tmp'});
    uabs=load(char(fil));
    uriem=load(char(fil2));
    ut=[ut;uabs(1,field)];
    uriemt=[uriemt;uriem(1,field)];
  end
  udiff=abs(ut-uriemt);
  figure(1);
  [ah]=axes;
  semilogy(dt*(1:300)',udiff,'r-*','MarkerSize',markersize_s,'Linewidth',linewidth_s);
  grid;
  set(ah,'Fontsize', fontsize_s);
  set(ah,'XMinorGrid','off','YMinorGrid','off');
  xlabel('simulation time [secs]');
  ylabel('abs(u_{riemann} - u_{ulsar})   at left boundary  (x=0 m)');
%%  title('differences ULSAR -- Riemann based absorbent B.C.s for rectang. channel shape');
 %%  anot=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
 %%  oh=text(1,1.035,anot);
 %%  set(oh(1),'FontSize',fontsize_s);
 %%  if (field==2)
 %%    axis([0 Lx 0.98 1.45]);
 %%  else
 %%    axis([0 Lx -0.5 0.5]);
 %%  end
    
    fig=char(strcat({'YUV/riem-comp-diff2.'},int2str(i),{'.tiff'}));
    print('-dtiff',fig,'-zbuffer');
    fig2=char(strcat({'YUV/riem-comp-diff2.'},int2str(i),{'.jpg'}));
    print('-djpeg100',fig2,'-zbuffer');
    
 %%  pause;
  end

field=2;
duriemdxt=[];duabsdxt=[];
if (1)
  for i=0:5:299
    i
    fil=strcat({'STEPS-rect-nouref/canal1d.'},int2str(i),{'.tmp'});
    fil2=strcat({'STEPS-rect-riemann/canal1driem.'},int2str(i),{'.tmp'});
    uabs=load(char(fil));
    duabsdx=diff(uabs(1:nnod,field));
    pp=sqrt(sum(duabsdx.^2));
    duabsdxt=[duabsdxt;pp];


    uriem=load(char(fil2));
    duriemdx=diff(uriem(1:nnod,field));
    qq=sqrt(sum(duriemdx.^2));
    duriemdxt=[duriemdxt;qq];

  end
  figure(1);
  [ah]=axes;
  semilogy(dt*(1:5:300)',duabsdxt,'r--^','MarkerSize',markersize_s,'Linewidth',linewidth_s);
  hold on;
  semilogy(dt*(1:5:300)',duriemdxt,'b-.*','MarkerSize',markersize_s,'Linewidth',linewidth_s);
  grid;
%  axis([0 10 0. 0.13]);
  set(ah,'Fontsize', fontsize_s);
  set(ah,'XMinorGrid','off','YMinorGrid','off');
  xlabel('simulation time [secs]');
  ylabel('||dh/dx||_2');
  title('||dh/dx||_2 vs. time [secs]')
  ol=legend('ULSAR B.C.','Riemann Invariants based B.C.',1);
  set(ol,'FontSize',fontsize_s);
%%  title('differences ULSAR -- Riemann based absorbent B.C.s for rectang. channel shape');
 %%  anot=char(strcat({'real time: '},num2str(i*dt),{' secs. '}));
 %%  oh=text(1,1.035,anot);
 %%  set(oh(1),'FontSize',fontsize_s);
 %%  if (field==2)
 %%    axis([0 Lx 0.98 1.45]);
 %%  else
 %%    axis([0 Lx -0.5 0.5]);
 %%  end
    
    fig=char(strcat({'YUV/riem-comp-diff2.'},int2str(i),{'.tiff'}));
    print('-dtiff',fig,'-zbuffer');
    fig2=char(strcat({'YUV/riem-comp-diff2.'},int2str(i),{'.jpg'}));
    print('-djpeg100',fig2,'-zbuffer');
    
 %%  pause;
  end
