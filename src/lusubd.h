// -*- mode: C++ -*- 
// $Id: lusubd.h,v 1.3 2001/07/13 03:48:23 mstorti Exp $
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
public:
  PFMatLU(Dofmap *dofmap_,Darray *da);
};

#endif
