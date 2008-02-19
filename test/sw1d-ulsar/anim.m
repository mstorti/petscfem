%clear all;
nod=load('canal1d.nod.tmp');
nnod=size(nod,1)-2;
x=nod(1:nnod,1);
y=nod(1:nnod,2);
h_i=nod(1:nnod,3);
con=load('canal1d.con.tmp');
%con=[con ones(size(con,1),1)];
dt=0.005;
field=2;Lx=2;
for i=0:1000
  i
  fil=strcat({'STEPS/canal1d.'},int2str(i),{'.tmp'});
  u=load(char(fil));
  plot(x(1:nnod),u(1:nnod,field),'r*');
  if (field==2)
    axis([0 Lx 0.999 1.025]);
  else
    axis([0 Lx -0.5 0.5]);
  end
  pause;
end
