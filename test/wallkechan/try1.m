##__INSERT_LICENSE__
## $Id: try1.m,v 1.2 2003/01/08 15:49:05 mstorti Exp $
d=10;                           # Supplemental value added to the diagonal
                                # in order to obtain a non-singular matrix.
N=100;                          # Size of the reduced problem
m=10;                           # Number of eqs. to be modified
lascal=1;                       # 
ladiag=0;

A=(2*rand(N)-1)+d*eye(N);
b=rand(N,1);

nw=N+1-(m:-1:1);
Anew=(2*rand(N)-1)+d*eye(N);
Anew(1:N-m,:) = A(1:N-m,:);
bnew = rand(N,1);
bnew(1:N-m)=b(1:N-m);
xex = Anew\bnew;

AA=[A, [zeros(N-m,m);lascal*eye(m)];
    Anew(nw,:), ladiag*eye(m)];

x=zeros(N,1);
while 1
  r=-Anew*x+bnew;
  RR=zeros(N+m,1);
  RR(1:N-m)=r(1:N-m);
  RR(N+(1:m))=r(nw);
  xx=AA\RR;
  x = x + RR([(1:N-m)'; N+(1:m)']);
  err=merr(x-xex)
  if err<1e-12
    break;
  endif
endwhile

return 
global linear_atv_data

linear_atv_data.A = Anew;
[x, error, total_iters] = gmres(0*bnew, bnew,"latv",[1e-12 100]);
error=error';
linear_atv_data.A = AA;
[xl, errorl, total_itersl] = gmres(0*BB, BB,"latv",[1e-12 100]);
errorl=errorl';
