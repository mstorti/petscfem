// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
// $Id: absolay.h,v 1.1 2006/04/17 03:08:14 mstorti Exp $
#ifndef PETSCFEM_ABSOLAY_H
#define PETSCFEM_ABSOLAY_H

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
  AbsorbingLayer(AbsorbingLayerFF *adv_diff_ff_=NULL) :
    adv_diff_ff(adv_diff_ff_), volume_flag(0) {};
  /** Destructor. */
  ~AbsorbingLayer() {delete adv_diff_ff;}
  
  /// The assemble function for the elemset. 
  NewAssembleFunction new_assemble;
  /// The ask function for the elemset. 
  ASK_FUNCTION;
};

#endif
