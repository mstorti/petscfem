// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: gtemplates.h,v 1.4 2004/12/05 15:38:37 mstorti Exp $
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int 
Tetra2TetraSplitter_v[] = {4,7,6,0,
			   4,5,8,1,
			   5,6,9,2,
			   8,9,7,3,
			   4,6,7,8,
			   4,5,6,8,
			   5,9,6,8,
			   9,7,8,6,GeomObject::NULL_NODE};

class Tetra2TetraSplitterClass : public Splitter {
private:
  dvector<int> subobj_conn;
public:
  Tetra2TetraSplitterClass() {
    subobj_conn.a_resize(2,8,4);
    subobj_conn.set(Tetra2TetraSplitter_v);
  }
  int size(GeomObject::Type t) const { 
    if (t==GeomObject::TetraT) return 8;
    else return 0;
  };
  const int* nodes(GeomObject::Type t,int j) { 
    if (t==GeomObject::TetraT) 
      return &subobj_conn.e(j,0);
    else assert(0);
  }
};

extern Tetra2TetraSplitterClass 
Tetra2TetraSplitter;

#endif
