// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id$
#ifndef PETSCFEM_FLOW_REVERSAL_H
#define PETSCFEM_FLOW_REVERSAL_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Switches boundary conditions from inlet to outlet
    for flow reversal. */
class flow_reversal : public GenLoadNS { 
public: 
  flow_reversal() {}
  void start_chunk() {}
  void end_chunk() {}
  void element_hook(int element) {}
  void q(FastMat2 &u,FastMat2 &flux,FastMat2 &jac);
};

#endif
