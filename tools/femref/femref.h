// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.17 2004/11/22 12:50:21 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/linkgraph.h>
#include <src/generror.h>

class GeomObject {
public:
  /// Types of objects
  enum Type { NULL_TYPE=0, 
	      OrientedEdgeT, EdgeT, 
	      OrientedTriT, TriT,
	      OrientedTetraT, TetraT};
  static int NULL_NODE;
protected:
  class Template {
    friend class GeomObject;
    dvector<int> perms_v;
    int size_m,dim_m,nperms_m;
    GeomObject::Type type;
    const char *label;
    Template(int sz,GeomObject::Type t,
	     int dim_a,int nperms_a,
	     const int *perms,const char *label_a);
  };
  static Template EdgeTemplate, TriTemplate, 
    OrientedTetraTemplate;
  Template *go_template;
private:
  int cs, canonical;
  dvector<int> nodes;
public:
  int size() { return go_template->size_m; }
  int csum() { return cs; }
  int dim() { return go_template->dim_m; }
  int nperms() { return go_template->nperms_m; }
  void clear() { go_template=NULL; nodes.clear(); canonical=0; cs=0; }
  void init(Type t,const int *nodes_a);
  const int* perm(int perm_indx) {
    return &(go_template->perms_v.e(perm_indx,0)); 
  }
  Type type() { return go_template->type;}
  GeomObject() : canonical(0), go_template(NULL) { }
  GeomObject(Type t,const int *nodes_a);
  void make_canonical();
  /// Compare two objects
  bool equal(GeomObject &go);
  void print(const char*s = NULL);
};

#endif
