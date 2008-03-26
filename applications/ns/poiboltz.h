// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharm.h,v 1.2 2002/12/16 04:11:35 mstorti Exp $
#ifndef POIBOLTZ_H
#define POIBOLTZ_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class poisson_boltzmann : public adaptor_pg {

  /* 
   *    Laplacian(phi) = A*sinh(B*phi)
   *
   *    A = 2*(ninf*z*F)/(eps*eps0)
   *    B = (z*F)/(R*Tabs)
   *
   */

  //o bulk concentration
  double ninf;
  //o valence
  int    z;
  //o valence
  //o relative permittivity
  double eps;
  //o vacuum permittivity
  double eps0;
  //o Faraday constant
  double F;
  //o absolute temperature
  double Tabs;
  //o ideal gas constant
  double R;

  // EDL tickness (if given, all values above are ignored)
  double Debye_length;

  //o A = 2*(ninf*z*F)/(eps*eps0)  [computed at elemset_init()]
  double A;
  //o B = (z*F)/(R*Tabs)           [computed at elemset_init()]
  double B;

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
