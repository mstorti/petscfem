// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: interplns.h,v 1.1 2003/06/09 02:37:15 mstorti Exp $
#ifndef PETSCFEM_INTERPLNS_H
#define PETSCFEM_INTERPLNS_H

#include <src/interpola.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Interpolates the elemset values for a given set of points. 
class interpolation_ns : public interpolation { 
public:
  /// The ask function for the elemset. 
  int ask(const char *jobinfo,int &skip_elemset) {
    skip_elemset = 1;
    DONT_SKIP_JOBINFO(comp_mat);
    DONT_SKIP_JOBINFO(comp_res);
    DONT_SKIP_JOBINFO(comp_mat_res);
    return 0;
  }
};

#endif
