// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: gtemplates.h,v 1.9 2005/01/03 03:15:22 mstorti Exp $
#ifndef PETSCFEM_GTEMPLATES_H
#define PETSCFEM_GTEMPLATES_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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
    if (t==GeomObject::OrientedTetraT) return 1;
    if (t==GeomObject::OrientedTriT) return 4;
    else return 0; 
  }
  const int* nodes(GeomObject::Type t,int j) const { 
    if (t==GeomObject::OrientedTriT) 
      return faces+3*j;
    else assert(0);
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class EdgeRefNodeTemplateClass 
  : public GeomObject::Template {
private:
  static int perm_v[];
public:
  EdgeRefNodeTemplateClass() 
    : GeomObject::Template(2,GeomObject::EdgeRefNodeT, 
			   1,1,perm_v,
			   "EdgeRefNode") { }
  int size(GeomObject::Type t) const { 
    return 0; 
  }
  const int* 
  nodes(GeomObject::Type t,int j) const { 
    assert(0);
  }
};

extern EdgeRefNodeTemplateClass 
EdgeRefNodeTemplate;

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

class Tetra2TetraSplitterClass : public Splitter {
private:
  dvector<int> subobj_conn;
  dvector<int> ref_nodes;
public:
  Tetra2TetraSplitterClass();
  int size(GeomObject::Type t) const;
  const int* nodes(GeomObject::Type t,int j) const;
  int size() const;
  const int* nodes(int j,GeomObject::Type &t) const;
  int nref_nodes() const;
  void 
  ref_node(int indx,
	   const GeomObject::Template *&tmpl,
	   int &nnod, const int *&nodes) const;
};

extern Tetra2TetraSplitterClass 
Tetra2TetraSplitter;

#endif
