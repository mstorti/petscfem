clear;
source("igrec.data.m.tmp");

x=(-L:ds:0.0);
p1=length(x);
xnod = [x' x' 0.*ones(p1,1)];
x=(-L:ds:0.0-ds);
p2=length(x)+p1;
p=p2-p1;
xnod = [xnod;x' -1.0.*x' 0.*ones(p,1)];
x=(ds:ds:L);
p3=length(x)+p2;
p=p3-p2;
xnod = [xnod; x' 0.*ones(p,1) 0.*ones(p,1)];
nnod = size(xnod,1);

fid=fopen("igrec.nod.tmp","w");
for j=1:nnod
  fprintf(fid,"%f %f %f\n",xnod(j,1),xnod(j,2),xnod(j,3));
endfor
fprintf(fid,"%f %f %f\n",-100-ds,-100-ds,0.);
fprintf(fid,"%f %f %f\n",-100-ds,+100+ds,0.);
fprintf(fid,"%f %f %f\n",+100+ds,0.,0.);
fprintf(fid,"%f %f %f\n",-100-ds*.5,-100-ds*.5,0.);
fprintf(fid,"%f %f %f\n",-100-ds*.5,+100+ds*.5,0.);
fprintf(fid,"%f %f %f\n",+100+ds*.5,0.,0.);
fclose(fid);

icone1=[1:p1-1]';
icone2=[2:p1]';
icone=[icone1 icone2];
icone1=[];icone2=[];
icone1=[p1+1:p2-1]';
icone2=[p1+2:p2]';
icone=[icone; icone1 icone2; p2 p1];
icone = [icone; nnod+4 1; nnod+5 p1+1];
p=size(icone,1);
fid=fopen("igrec1.con.tmp","w");
for k=1:p
  fprintf(fid,"%d %d\n",icone(k,1),icone(k,2));
endfor
fclose(fid);

icone1=[];icone2=[];icone=[];
icone1=[p2+1:p3-1]';
icone2=[p2+2:p3]';
icone=[p1 p2+1;icone1 icone2; nnod nnod+6];
p=size(icone,1);
fid=fopen("igrec2.con.tmp","w");
for k=1:p
  fprintf(fid,"%d %d\n",icone(k,1),icone(k,2));
endfor
fclose(fid);

h_ini=h0*ones(nnod+6,1);
h_ini(1)=h0*1.5;
h_ini(p1+1)=h0*2.;
h_ini(nnod+1)=0.;
h_ini(nnod+2)=0.;
h_ini(nnod+3)=0.;
u_ini=u0*ones(nnod+6,1);
u_ini(nnod+1)=0.;
u_ini(nnod+2)=0.;
u_ini(nnod+3)=0.;
fid=fopen("igrec.ini.tmp","w");
for i=1:nnod+6
fprintf(fid,"%f %f\n",u_ini(i),h_ini(i));
endfor
fclose(fid);

if (1)
rest = [1 2 3 4 nnod+1 nnod+4 -1.];
rest = [rest; p1+1 p1+2 p1+3 p1+4 nnod+2 nnod+5 -1.];
fid=fopen("igrec1.rest.tmp","w");
fprintf(fid,"%d %d %d %d %d %d %f\n",rest');
fclose(fid);
rest=[];
rest=[nnod nnod-1 nnod-2 nnod-3 nnod+3 nnod+6  1.];
fid=fopen("igrec2.rest.tmp","w");
fprintf(fid,"%d %d %d %d %d %d  %f\n",rest');
fclose(fid);
endif

if (1)
fixa=[];
fixa = [fixa; nnod+4 1 u0+0.*u0; nnod+4 2 h0+0.*h0;nnod+5 1 u0+0.*u0;nnod+5 2 h0+0.*h0;nnod+6 1 u0;nnod+6 2 h0];
fid=fopen("igrec.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);
endif

clear all;