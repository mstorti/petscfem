function err = checksym(vec)

  estela = 79*30+(1:20);
  for k=1:4
    rho=vec(:,k);
    rho(2449:-1:2430) = rho(estela);
    rho=reshape(rho,79,31);
    sig=1;
    if k==3 
      sig=-1;
    end
    err(k) = max(max(abs(rho-sig*rho(79:-1:1,:))));
  end
endfunction
