## $Id: proc2.m,v 1.1 2003/09/18 14:32:52 mstorti Exp $
hhh=readconv("Residual_norm","nohup.out");
hh=readconv("norm_res","nohup.out");
h=readconv("delta_u","nohup.out");

np = 3000;
n2=length(hhh);
n1=max([n2-np+1,1]);

printf("%d read, %d plotted\n",length(hhh),np);

semilogy(hhh(n1:n2))
pause

nnwt=2;
nstep = fix(length(hh)/2);
hh=reshape(hh(1:2*nstep),2,nstep)';
semilogy(hh)
pause

plot(h)


