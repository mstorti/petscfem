// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: invcoupl.h,v 1.3 2003/02/27 03:32:41 mstorti Exp $
#ifndef PETSCFEM_INVCOUPL_H
#define PETSCFEM_INVCOUPL_H

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
