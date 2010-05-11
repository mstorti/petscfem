// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharm.h,v 1.2 2002/12/16 04:11:35 mstorti Exp $
#ifndef charge_cons_H
#define charge_cons_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class charge_cons : public adaptor_pg {
private:

  string velname;
  double* velptr;
  FastMat2 velcol;
  FastMat2 vel;

  string sigma_name;
  double* sigmaptr;
  FastMat2 sigmacol;
  FastMat2 sigma;

  string gamma_name;
  double* gammaptr;
  FastMat2 gammacol;
  FastMat2 gamma;

  string beta_name;
  double* betaptr;
  FastMat2 betacol;
  FastMat2 grad_beta;
  FastMat2Tmp tmp;

public:
  void elemset_init();
  void elem_init();
  void pg_connector(const FastMat2 &xpg,
		    const FastMat2 &state_old_pg,
		    const FastMat2 &grad_state_old_pg,
		    const FastMat2 &state_new_pg,
		    const FastMat2 &grad_state_new_pg,
		    FastMat2 &res_pg,FastMat2 &mat_pg);

public:
  Nodedata *nodedata;
  ASSEMBLE_FUNCTION;
};

#endif
