## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: procfs.m,v 1.2 2002/07/17 02:55:01 mstorti Exp $

## Author: Mario Storti
## Keywords: fractional step, square cavity
## usage: 
source("data.m.tmp");
xnod = aload("sqcav.nod.tmp");
y=xnod(1:N+1,2);
x=xnod(1:N+1:(N+1)^2,1);

U=aload("sqcav.weak_form_1.tmp");
u=reshape(U(:,1),N+1,N+1);
v=reshape(U(:,2),N+1,N+1);
p=reshape(U(:,3),N+1,N+1);
