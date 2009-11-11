// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharm.h,v 1.2 2002/12/16 04:11:35 mstorti Exp $
#ifndef POISSON_H
#define POISSON_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class poisson: public adaptor_pg {

  
  double scale_factor;

  //o relative permittivity
  double eps;
  //o vacuum permittivity
  double eps0;
  //o Faraday constant
  double F;
  double A;
  
  string concentration_name;
  double* concptr;
  FastMat2 concol;
  FastMat2 tmp;
  FastMat2 conc;

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
