// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nsitetlesf.h,v 1.1.2.1 2007/02/27 21:10:32 mstorti Exp $
#ifndef PETSCFEM_NSITETLESF_H
#define PETSCFEM_NSITETLESF_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_full : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

#endif
