// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: invcoupl.h,v 1.2 2003/02/24 00:14:23 mstorti Exp $
#ifndef PETSCFEM_INVCOUPL_H
#define PETSCFEM_INVCOUPL_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class FastMat2Tmp {
private:
  vector<FastMat2 *> store_v;
  FastMat2 **store;
  int size;
  void sync() {
    size = store_v.size();
    store = store_v.begin();
  }
public:
  FastMat2Tmp() : size(0), store(NULL) {}
  FastMat2 & operator()(int j) {
    if (j>=size) { 
      store_v.resize(j+1); 
      for (int k=size; k<j; k++) store_v[k] = NULL;
      sync();
    }
    if (!store[j]) store[j] = new FastMat2;
    return *store[j];
  }
  void clear() {
    for (int j=0; j<size; j++) delete store[j];
    store_v.clear();
    sync();
  }
  ~FastMat2Tmp() { clear(); }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class inviscid_coupling : public adaptor_pg {
  double viscosity;
  FastMat2Tmp tmp;
  // FastMat2 tmp1,tmp2,tmp3,tmp4,gsopg;
public:
  void elemset_init();
  void elemset_end();
  void pg_connector(const FastMat2 &xpg,
		    const FastMat2 &state_old_pg,
		    const FastMat2 &grad_state_old_pg,
		    const FastMat2 &state_new_pg,
		    const FastMat2 &grad_state_new_pg,
		    FastMat2 &res_pg,FastMat2 &mat_pg);
};

#endif
