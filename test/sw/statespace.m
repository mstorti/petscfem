clear;
wm=(-4.5:0.1:4.5)';wp=(-4.5:0.1:4.5)';
h=(0.8:0.01:1)';

for i=1:length(wm)
  um(:,i)=wm(i)+4*sqrt(h);
  up(:,i)=wp(i)-4*sqrt(h);
endfor
um=um';up=up';
hold off;
figure(1);
axis([0.9 1 0.1 0.5])
plot(h,up);hold on;plot(h,um);
hold off