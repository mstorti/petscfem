// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: interpola.h,v 1.1 2003/06/08 21:42:01 mstorti Exp $
#ifndef PETSCFEM_INTERPOLA_H
#define PETSCFEM_INTERPOLA_H

#include <src/pfobject.h>
#include <src/dvector.h>

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class interpolation_e : public NewElemset { 
public:

  interpolation_e() {}
  ~interpolation_e() {}
  /// The ask function for the elemset. 

  int ask(const char *jobinfo,int &skip_elemset);

  virtual void initialize();

  void new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
		    const Dofmap *dofmap,const char *jobinfo,
		    const ElementList &elemlist,
		    const TimeData *time_data);

};

#endif
