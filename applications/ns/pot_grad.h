// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: qharm.h,v 1.2 2002/12/16 04:11:35 mstorti Exp $
#ifndef pot_grad_H
#define pot_grad_H

#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class pot_grad : public adaptor_pg {
private:
 
  // FastMat2 Dm;
  //FastMat2 Cm;
 



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
