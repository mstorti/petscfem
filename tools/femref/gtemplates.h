// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: gtemplates.h,v 1.3 2004/11/23 18:41:10 mstorti Exp $
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
			   3,11,perm_v,
			   "OrientedTetra") { }
  int size(GeomObject::Type t) const { 
    if (t==GeomObject::OrientedTriT) return 4;
    else return 0; 
  }
  const int* nodes(GeomObject::Type t,int j) const { 
    assert(t==GeomObject::OrientedTriT); 
    return faces+3*j;
  }
};

extern OrientedTetraTemplateClass 
OrientedTetraTemplate;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class OrientedTriTemplateClass 
  : public GeomObject::Template {
private:
  static int perm_v[];
public:
  OrientedTriTemplateClass() 
    : GeomObject::Template(3,GeomObject::OrientedTriT, 
			   2,5,perm_v,
			   "OrientedTri") { }
  int size(GeomObject::Type t) const { }
  const int* nodes(GeomObject::Type t,int j) const { }
};

extern OrientedTriTemplateClass 
OrientedTriTemplate;

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define TEMPLATE(TYPE,size,dim,nperms)				\
GeomObject::Template						\
GeomObject::TYPE##Template(size,GeomObject::TYPE##T,		\
			 dim,nperms,TYPE##TemplatePerm_v,#TYPE)

static int EdgeTemplatePerm_v[] = {1,0,GeomObject::NULL_NODE};
TEMPLATE(Edge,2,1,1);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static int 
TriTemplatePerm_v[] = {1,2,0,
		       2,0,1,
		       0,2,1,
		       2,1,0,
		       1,0,2,GeomObject::NULL_NODE};
TEMPLATE(Tri,3,2,5);
#endif

#endif