function w=zhomo2(b,nx,ny,str)
%function w=zhomo2(b,nx,ny,str)

if nargin==3
 str=[1 0 1   1 0 1];
end
strx=str(1:3);
stry=str(4:6);

wx=b(1)+(onedstr(strx,nx))'*(b(2)-b(1));
wy=b(3)+onedstr(stry,ny)*(b(4)-b(3));

keyboard;pause

w=kron(wx,ones(size(wy)))+i*kron(ones(size(wx)),wy);
