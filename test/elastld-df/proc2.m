###key proc2.m
### $Id: $
source("data.m.tmp");

q=aload("elastld-some-rslt.tmp"); 
nt = rows(q);
# plot((1:nt)'*Dt,q(:,3));
plot((1:nt)'*Dt,q(:,2:4));
