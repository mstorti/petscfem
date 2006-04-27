clear all;

M=21;
N=21;

load xnod.out;
load icone.out;
load state.out;
icone = icone+1;

u = state(:,1);
v = state(:,2);
p = state(:,3);

x = linspace(0,1,M);
y = linspace(0,1,N);
[X,Y] = meshgrid(x,y);
U = reshape(u,M,N)';
V = reshape(v,M,N)';
P = reshape(p,M,N)';


streamline(X,Y,U,V, [x(1:2:end),1-x(1:2:end)], [y(1:2:end),y(1:2:end)]);
axis equal; grid on;

%quiver(X,Y,U,V, 0);
%axis equal;

%surf(X,Y,P)
%axis equal;
