// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharm.h,v 1.2 2002/12/16 04:11:35 mstorti Exp $
#ifndef electrophoresisM2_H
#define electrophoresisM2_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class electrophoresisM2 : public adaptor_pg {
private:
  FastMat2 diff;
  FastMat2 zeff;
  FastMat2 adveff;
  FastMat2 K;
  FastMat2 K1;
  FastMat2 Dm;
  FastMat2 Cm;
  FastMat2 Hm;
  FastMat2 ndof_aux;
  FastMat2 A_supg;
  FastMat2 B_supg;
  FastMat2 C_supg;
  FastMat2 D_supg;
  FastMat2 tau_supg;
  double mu;
  double supg_fact;
  double r_fact;
  double m_fact;

  string velname;
  double* velptr;
  FastMat2 velcol;
  FastMat2 vel;

  string potname;
  double* potptr;
  FastMat2 potcol;
  FastMat2 pot;
  FastMat2 celec;
  FastMat2 lapfi;
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
