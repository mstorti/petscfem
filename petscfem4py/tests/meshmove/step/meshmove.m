clear all;
load xnod.out;
load icone.out;
icone = icone+1;
trimesh(icone, xnod(:,1), xnod(:,2))
