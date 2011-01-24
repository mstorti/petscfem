## $Id: mkmesh.m,v 1.3 2006/03/12 12:07:37 mstorti Exp $
source("data.m.tmp");

w = zhomo([0,1,0,1],Nr+1,Nphi+1);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1,4,3,2]);
r = R+(2*xnod(:,1)-1)*DR/2;
phi = xnod(:,2)*pi/2;

xnod = [r.*cos(phi),r.*sin(phi)];

asave("bubble.nod.tmp",xnod);
asave("bubble.con.tmp",icone);

extlay=(Nphi+1)*Nr+(1:(Nphi))';
asave("bubble.load.tmp",[extlay,extlay+1]);

tol = 1e-6;
tmp = find(abs(phi)<tol);
pffixa("bubble.fixa.tmp",tmp,2,0.0);

tmp = find(abs(phi-pi/2)<tol);
pffixa("bubble.fixa.tmp",tmp,1,0.0,"a");
