// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: gtemplates.h,v 1.1 2004/11/22 16:07:53 mstorti Exp $
#ifndef PETSCFEM_GTEMPLATES_H
#define PETSCFEM_GTEMPLATES_H

class OrientedTetraTemplateClass 
  : public GeomObject::Template {
private:
  static int faces[];
  static int perm_v[];
public:
  OrientedTetraTemplateClass() 
    : GeomObject::Template(4,GeomObject::OrientedTetraT, 
			   3,11,OrientedTetraTemplatePerm_v,
			   "OrientedTetra") { }
  int size(GeomObject::Type t) { 
    if (t==GeomObject::OrientedTriT) return 4;
    else return 0; 
  }
  const int* nodes(GeomObject::Type t,int j) { 
    assert(t==GeomObject::OrientedTriT); 
    return faces+3*j;
  }
};

GeomObject::Template *OrientedTetraTemplate;

#endif
