## $Id: proc.m,v 1.1 2005/02/25 01:59:47 mstorti Exp $

asave("mesh1.nod",[rand(4,2),zeros(4,1)]);
system("project.bin");
q = aload("qq.dat");

xeh = q(:,1:3);
xprojh = q(:,4:6);

plot(xeh(:,1),xeh(:,2),'og',xprojh(:,1),xprojh(:,2),'or',\
     [xeh(:,1)';xprojh(:,1)'], [xeh(:,2)';xprojh(:,2)'],'b');
