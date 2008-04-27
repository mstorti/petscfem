// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id$
#ifndef electrophoresis_H
#define electrophoresis_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class electrophoresis : public adaptor_pg {

  double diff;
  int    z;
  double el;
  double mu;
  double rm;
  double r;
  double N;

  string velname;
  double* velptr;
  FastMat2 velcol;
  FastMat2 vel;

  string potname;
  double* potptr;
  FastMat2 potcol;
  FastMat2 pot;
  FastMat2 celec;
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
