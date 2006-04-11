// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: errestim.h,v 1.2 2006/04/11 00:48:10 mstorti Exp $
#ifndef ERRESTIM_H
#define ERRESTIM_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class error_estimator : public adaptor_pg {
  double norm_expo;
  FastMat2 du, tmp;
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
