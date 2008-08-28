//__INSERT_LICENSE__
//$Id$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./nsi_tet.h"
#include "./nsitetlesfbf.h"


#define STANDARD_UPWIND
#define USE_FASTMAT

extern TextHashTable *GLOBAL_OPTIONS;

#define STOP {PetscFinalize(); exit(0);}
   
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

void nsi_tet_les_full_bf::
bf_init(Nodedata* nodedata)
{ 
  int ierr;

  TGETOPTDEF_S_ND(thash,string,potential_field,potential);
  TGETOPTDEF_S_ND(thash,string,charge_field,charge);

  bool found; int ncols;
  found = nodedata->get_field(potential_field,&ncols,&potn_ptr);
  assert(found); assert(ncols==1);
  found = nodedata->get_field(charge_field,&ncols,&chrg_ptr);
  assert(found); assert(ncols==1);

  potncol.resize(1,nel);
  chrgcol.resize(1,nel);
  grad_potn.resize(1,nodedata->ndim);

}

void nsi_tet_les_full_bf::
bf_eval_el(int iele)
{ 
  for (int i=0; i<nel; i++) {
    int     node = icone[iele*nel+i];
    double  potn = potn_ptr[node-1];
    double  chrg = chrg_ptr[node-1];
    potncol.setel(potn, i+1);
    chrgcol.setel(chrg, i+1);
  }
}

void nsi_tet_les_full_bf::
bf_eval_pg(FastMat2& SHAPE,FastMat2& DSHAPE, FastMat2& BODY_FORCE)
{
  // compute potential, charge and grad of potential
  double potn = tmp(1).prod(SHAPE,potncol,-1,-1);
  double chrg = tmp(1).prod(SHAPE,chrgcol,-1,-1);
  grad_potn.prod(DSHAPE,potncol,1,-1,-1);
  
  // compute body force
  BODY_FORCE
    .set(grad_potn)
    .scale(-chrg);
    ;

}
