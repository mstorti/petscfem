// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
// $Id: absolay.h,v 1.2 2006/04/18 02:16:56 mstorti Exp $
#ifndef PETSCFEM_ABSOLAY_H
#define PETSCFEM_ABSOLAY_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/elemset.h>
#include "./advective.h"

/** Perfectly Matched Layer type elemset */
class AbsorbingLayer : public NewElemset { 
protected:
  /** A pointer to the flux function. Different
      physics are obtained by redefining the flux function
  */
  NewAdvDifFF *adv_diff_ff;

public:
  /// Contructor from the pointer to the fux function
  AbsorbingLayer(NewAdvDifFF *adv_diff_ff_a=NULL) :
    adv_diff_ff(adv_diff_ff_a) {};
  /** Destructor. */
  ~AbsorbingLayer() {delete adv_diff_ff;}
  
  /// The assemble function for the elemset. 
  NewAssembleFunction new_assemble;
  /// The ask function for the elemset. 
  ASK_FUNCTION;
};

#endif
