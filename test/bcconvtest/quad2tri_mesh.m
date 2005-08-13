function icone=quad2tri_mesh(icone)
%
%   quad2tri_mesh(icone)
%

[numel,nen]=size(icone);

iconew=[];
for k=1:numel,
  vvv  = icone(k,:);
  ico  = q2t(vvv);
  iconew = [iconew;ico(1,:);ico(2,:)];
end
icone=iconew;
