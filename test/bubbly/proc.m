if 0
  petsc_data_name="a_ll_000";
  petscload
  A=getblock(a_ll_000);
else
  petsc_data_name="mat.output";
  petscload
  A=getblock(A);
  AA=getblock(AA);
endif
