// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharmm.h,v 1.2 2002/12/17 02:23:40 mstorti Exp $
#ifndef QHARMM_H
#define QHARMM_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Quasi-harmonic equation, for multiple degrees of freedom
class qharmm : public adaptor_pg {
private:
    FastMat2 cond,tmp1,tmp2;
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
