// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: interpola.h,v 1.2 2003/06/08 22:51:03 mstorti Exp $
#ifndef PETSCFEM_INTERPOLA_H
#define PETSCFEM_INTERPOLA_H

//#include <src/pfobject.h>
#include <src/dvector.h>

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This objects provides contranints that restrict the velocity
    field at a certain set of nodes to be irrotational.  */ 
class interpolation : public BasicObject {
private:
  NewElemset *volume_elemset;
  Nodedata *nodedata;
  dvector<double> points;
  int npoints;
public:
  /** Ctor.  */ 
  interpolation() {}
  /** Dtor.  */ 
  ~interpolation() {}
  /** Specific read function for this obkect.  */ 
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap);
};
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Interpolates the elemset values for a given set of points. 
class interpolation : public NewElemset { 
private:
  /// The 
  Nodedata *nodedata;
  dvector<double> points;
  int npoints;
public:
  interpolation() {}
  ~interpolation() {}
  /// The ask function for the elemset. 
  int ask(const char *jobinfo,int &skip_elemset);
  virtual void initialize();
  void new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
		    const Dofmap *dofmap,const char *jobinfo,
		    const ElementList &elemlist,
		    const TimeData *time_data);
  void read(FileStack *fstack);
  /// 
  virtual char *jobinfo_stage();
};

#endif
