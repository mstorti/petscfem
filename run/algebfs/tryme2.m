output
d=max(abs(D));
Iu = find(d==0);
Ip = find(d>0);

clear d
disp('A:'); A=takef(A,Iu,Iu);
disp('B:'); B=takef(B,Iu,Ip);
disp('C:'); C=takef(C,Ip,Iu);
disp('D:'); D=takef(D,Ip,Ip);

disp('resp:'); resp=takef(resp,Ip);
disp('resu:'); resu=takef(resu,Iu);
disp('resH:'); resH=takef(resH,Ip);
disp('dp:'); dp=takef(dp,Ip);
disp('du:'); du=takef(du,Iu);

dxx=[A B; C D]\[resu;resp];
duu=dxx(1:722);
dpp=dxx(722+(1:440));
dx=[du;dp];
