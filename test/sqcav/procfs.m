## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: procfs.m,v 1.6 2004/10/25 02:43:27 mstorti Exp $

## Author: Mario Storti
## Keywords: fractional step, square cavity
## usage: 
source("data.m.tmp");
ghia
xnod = aload("sqcav.nod.tmp");
y=xnod(1:N+1,2);
x=xnod(1:N+1:(N+1)^2,1);

U=aload("outvector0.out");
u=reshape(U(:,1),N+1,N+1);
v=reshape(U(:,2),N+1,N+1);
p=reshape(U(:,3),N+1,N+1);
ghia;

plot([u(:,N/2+1),ulu(:,N/2+1)],y,ug_1000,yg,'+')
