// -*- mode: C++ -*- 
// $Id: lusubd.h,v 1.2 2001/07/13 01:30:38 mstorti Exp $
#ifndef LUSUBD_H
#define LUSUBD_H

#include <vector>
#include "libretto.h"
#include "dofmap.h"

class PFMat {
public:
  virtual void set_profile(Darray *da)=0;
  virtual ~PFMat()=0;
};

class PFMatLU : public PFMat {
  Dofmap *dofmap;
public:
  PFMatLU(Dofmap *dofmap_,Darray *da);
};

#endif
