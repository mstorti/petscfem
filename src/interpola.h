// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: interpola.h,v 1.3 2003/06/09 02:37:18 mstorti Exp $
#ifndef PETSCFEM_INTERPOLA_H
#define PETSCFEM_INTERPOLA_H

#include <src/dvector.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Interpolates the elemset values for a given set of points. 
class interpolation : public NewElemset { 
private:
  dvector<double> points, interp_coef;
  dvector<int> point2elem;
  int npoints;
  int initialized;
  int ndim;
public:
  interpolation() : initialized(0) {}
  ~interpolation() {}
  // void initialize();
  void new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
		    const Dofmap *dofmap,const char *jobinfo,
		    const ElementList &elemlist,
		    const TimeData *time_data);
  void read(FileStack *fstack);
};

#endif
