// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.15 2004/11/21 21:22:42 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/linkgraph.h>
#include <src/generror.h>

#if 0
typedef unsigned int Node;

class GeomObjectBasic;
class GeomObjectId;

class GeomObject {
private:
  GeomObjectBasic *gobj;
public:
  /// Types of objects
  enum Type { NULL_TYPE=0, OrientedEdgeT, EdgeT, OrientedTriT, TriT };
  /// Ctor of empty object
  GeomObject();
  /// Set object to especific instance
  void init(int sz,const int *nodes,Type t);
  /// Constructor
  GeomObject(int sz,const int *nodes,
	     Type t) { init(sz,nodes,t); }
  /// Compare two objects
  bool equal(GeomObject &go);
  /// Dimension of the object
  int dim();
  /** Check sum of the object, may be used as 
      a weaker comparison operator. */ 
  int csum();
  /// Number of objects of type #t# adjacent to this object. 
  int size(Type t);
  /// Number of nodes adjacent to this object. 
  int size();
  /** Number of permutations that leave 
      this object invariant. */ 
  int nperms();
  /** Returns the #perm_indx# permutatation for 
      this kinf of object. */
  const int* perm(int perm_indx);
  /// Type of this object. 
  Type type();
  /// 
  void get_adjacency(const GeomObjectId &id1,
		     GeomObject::Type t,
		     std::list<GeomObjectId> &adjacent) { }
  /// Bring the node sequence to canonical order
  void make_canonical();
  /// Clears the object and set to the empty object. 
  void clear();
  /// Dtor. 
  ~GeomObject();
};

/** Geometric objects are stored via a kind of
    smart pointer. This is the real object stored 
    via an indirect pointer in a #GeomObject#. */ 
class GeomObjectBasic {
public:
  /// Compares two objects
  virtual bool equal(GeomObjectBasic &)=0;
  /// Returns the dimension of the object
  virtual int dim()=0;
  /// Returns the a fast check sum for the object. 
  virtual int csum()=0;
  /// Returns the number of nodes connected to this object. 
  virtual int size();
  /** Returns the number of object of type #t#
      connected to #*this#. */ 
  virtual 
  int size(GeomObject::Type t) {return 0; }
  /// Type of this object. 
  virtual GeomObject::Type type() { 
    return GeomObject::NULL_TYPE; 
  }
  /** Number of permutations that leave 
      this object invariant. */ 
  virtual int nperms() { return 0; }
  /** Returns the #perm_indx# permutatation for 
      this kinf of object. */
  virtual const int* perm(int perm_indx);
  /// Creates objects on demand. 
  static GeomObjectBasic 
  *factory(int sz,const int *nodes, GeomObject::Type t);
};
#endif

class GeomObjectSeq : public GeomObjectBasic {
protected:
  class Template {
    friend class GeomObjectSeq;
    static int NULL_NODE;
    dvector<int> perms_v;
    int size_m,dim_m,nperms_m;
    GeomObject::Type type;
    Template(int sz,GeomObject::Type t,
	       int dim_a,int nperms_a,const int *perms);
  };
  static Template EdgeTemplate;
  Template *go_template;
private:
  int cs, canonical;
  dvector<int> nodes;
public:
  int size() { return go_template->size_m; }
  int csum() { return cs; }
  int dim() { return go_template->dim_m; }
  int nperms() { return go_template->nperms_m; }
  const int* perm(int perm_indx) {
    return &(go_template->perms_v.e(perm_indx,0)); 
  }
  GeomObject::Type type() { return go_template->type; }
  GeomObjectSeq(GeomObject::Type t,const int *nodes_a);
  void make_canonical();
  /// Compare two objects
  bool equal(GeomObject &go);
};

#if 0
#define GO_CLASS(TYPE)					\
class TYPE : public GeomObjectSeq {			\
public:							\
  TYPE(const int *nodes_a)				\
    : GeomObjectSeq(&GeomObjectSeq::TYPE##Template,	\
		    nodes_a) { }			\
}

GO_CLASS(Edge);
#endif

#endif
