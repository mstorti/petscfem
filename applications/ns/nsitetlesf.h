// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nsitetlesf.h,v 1.1.4.1 2007/03/08 19:22:16 dalcinl Exp $
#ifndef PETSCFEM_NSITETLESF_H
#define PETSCFEM_NSITETLESF_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_full : public ns_volume_element { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

#endif
