if n_int_v(rows(n_int_v)) != rows(A)
  error("Sizes do not match!!");
endif

neq = n_int_v(rows(n_int_v));
A_iisd=zeros(neq);
for k=1:nproc
  indx = [(n_loc_v(k)+1):n_loc_v(k+1)];
  eval(sprintf("A_iisd(indx,indx)=a_ll_%03d;",k-1));
endfor

n_int=n_int_v(rows(n_int_v))-n_int_v(1);
n_loc=n_loc_v(rows(n_loc_v))-n_loc_v(1);
II = n_loc+(1:n_int)';
LL = (1:n_loc)';

A_iisd(LL,II) = A_LI;
A_iisd(II,LL) = A_IL;
A_iisd(II,II) = A_II;

A_iisd = A_iisd(map,:);
A_iisd = A_iisd(:,map);
