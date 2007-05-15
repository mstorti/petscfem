// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharm.h,v 1.2 2002/12/16 04:11:35 mstorti Exp $
#ifndef POIBOLTZ_H
#define POIBOLTZ_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class poisson_boltzmann : public adaptor_pg {

  double ninf;
  int    z;
  double eps;
  double eps0;
  double kb;
  double Tabs;

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
