fr=sqrt(U(:,2).^2+U(:,3).^2)./U(:,1)./sqrt(gravity*U(:,1));
fr=reshape(fr,Nx+1,Ny+1);
II=[1:Ny+1 Ny:-1:2 1:Ny+1 Ny:-1:2 1:Ny+1];
