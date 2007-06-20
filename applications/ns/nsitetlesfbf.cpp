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

  TGETOPTDEF_S_ND(thash,string,pot1name,pot1);
  TGETOPTDEF_S_ND(thash,string,pot2name,pot2);

  SGETOPTDEF_ND(double, param1, 0.0 );
  SGETOPTDEF_ND(double, param2, 0.0 );

  bool found; int ncols;
  found = nodedata->get_field(pot1name,&ncols,&field1);
  assert(found); assert(ncols==1);
  found = nodedata->get_field(pot2name,&ncols,&field2);
  assert(found); assert(ncols==1);

  pot1vec.resize(nel);
  pot2vec.resize(nel);

  pot1col.resize(1,nel);
  pot2col.resize(1,nel);

  grad_pot1.resize(1,nodedata->ndim);
  grad_pot2.resize(1,nodedata->ndim);

}

void nsi_tet_les_full_bf::
bf_eval_el(int iele)
{ 
  for (int i=0; i<nel; i++) {
    int node = icone[iele*nel+i];
    pot1vec[i] = field1[node-1];
    pot2vec[i] = field2[node-1];
  }
  pot1col.set(&pot1vec[0]);
  pot2col.set(&pot2vec[0]);
}

void nsi_tet_les_full_bf::
bf_eval_pg(FastMat2& SHAPE,FastMat2& DSHAPE, FastMat2& BODY_FORCE)
{
  // compute potentials
  double pot1 = tmp(1).prod(SHAPE,pot1col,-1,-1);
  double pot2 = tmp(1).prod(SHAPE,pot2col,-1,-1);

  // compute potential gradients
  grad_pot1.prod(DSHAPE,pot1col,1,-1,-1);
  grad_pot2.prod(DSHAPE,pot2col,1,-1,-1);
  
  // compute body force
  BODY_FORCE
    .axpy(grad_pot1, param1)
    .axpy(grad_pot2, param2)
    ;

}
