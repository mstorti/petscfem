n=20;

fid = fopen("matrix1.dat","w");
fprintf(fid,"%d\n",n-1);

h=1/n;
for j=0:n-2;
  fprintf(fid,"%g\n",h^2);
endfor

A = zeros(n-1);
for j=0:n-2
  if j!=0, 
    fprintf(fid,"%d %d %g\n",j,j-1,-1); 
    A(j+1,j) = -1;
  endif
  fprintf(fid,"%d %d %g\n",j,j,2); 
  A(j+1,j+1) = 2;
  if j!=n-2, 
    fprintf(fid,"%d %d %g\n",j,j+1,-1); 
    A(j+1,j+2) = -1;
  endif
endfor

max(eig(A))

fclose(fid);
