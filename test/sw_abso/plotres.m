load igrec0.out.tmp;
load igrec.nod.tmp;
x=igrec;
load igrec.ini.tmp;
ini=igrec;

h=reshape(igrec0(:,2),size(x,1),size(igrec0,1)/size(x,1));
u=reshape(igrec0(:,1),size(x,1),size(igrec0,1)/size(x,1));
clear igrec0;

M = moviein(size(h,1)+1);
plot3(x(1:size(x,1)-6,1),x(1:size(x,1)-6,2),ini(1:size(x,1)-6,2),'.');
%axis([-100 100 -100 100 0.99 1.1]);
M(:,1) = getframe;
for i=1:size(h,2),...
plot3(x(1:size(x,1)-6,1),x(1:size(x,1)-6,2),h(1:size(x,1)-6,i),'.');
%axis([-100 100 -100 100 0.99 1.1]);
M(:,i+1) = getframe;
end
movie(M,1);
clear all;