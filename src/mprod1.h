//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & FastMat2::prod(__MATS__
                          const int m,INT_VAR_ARGS_ND) {

  vector<const FastMat2 *> mat_list;
  __PACK__
  
  Indx indx;
  indx.push_back(m);
  READ_INT_ARG_LIST(indx);
  prod(mat_list,indx);
  return *this;
}
