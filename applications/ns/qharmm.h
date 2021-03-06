// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharmm.h,v 1.6 2006/10/23 02:43:18 mstorti Exp $
#ifndef PETSCFEM_QHARMM_H
#define PETSCFEM_QHARMM_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Quasi-harmonic equation, for multiple degrees of freedom
class qharmm : public adaptor_pg {
private:
  double Dt,rec_Dt;
  FastMat2 cond,C,Cp,x_ref,G;
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
