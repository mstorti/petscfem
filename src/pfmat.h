// -*- mode: C++ -*- 
// $Id: pfmat.h,v 1.1 2001/07/14 13:39:38 mstorti Exp $
#ifndef LUSUBD_H
#define LUSUBD_H

#include <vector>
#include "libretto.h"
#include "fem.h"
#include "dofmap.h"

class PFMat {
public:
  virtual ~PFMat()=0;
};

class PFMatLU : public PFMat {
  Dofmap *dofmap;
  int n_int,n_loc;
  vector<int> map;
  Mat A_LL,A_LI,A_IL,A_II;
public:
  PFMatLU(Dofmap *dofmap_,Darray *da);
};

#endif
