function [ico_t] = q2t(icone)
%
%  convert quadrangular mesh in triangular mesh
%
%         [ico_t] = q2t(icone)
%

[nele,nen] = size(icone);

if nen~=4, error(' q2t Error [000]');end
ico_t = icone(:,[1,2,3]);
ico_t = [ico_t ; icone(:,[1,3,4])];