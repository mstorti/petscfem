// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharmm.h,v 1.3 2003/03/22 22:20:32 mstorti Exp $
#ifndef PETSCFEM_QHARMM_H
#define PETSCFEM_QHARMM_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Quasi-harmonic equation, for multiple degrees of freedom
class qharmm : public adaptor_pg {
private:
  FastMat2 cond,C,x_ref;
  FastMat2Tmp tmp;
public:
  void elemset_init();
  void pg_connector(const FastMat2 &xpg,
		    const FastMat2 &state_old_pg,
		    const FastMat2 &grad_state_old_pg,
		    const FastMat2 &state_new_pg,
		    const FastMat2 &grad_state_new_pg,
		    FastMat2 &res_pg,FastMat2 &mat_pg);
};

#endif
