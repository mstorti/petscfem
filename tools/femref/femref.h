// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.13 2004/11/21 20:21:53 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/linkgraph.h>
#include <src/generror.h>

typedef unsigned int Node;

class GeomObjectBasic;

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

class GeomObjectSeq {
private:
  dvector<int> nodes;
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
  virtual Template *go_template()=0;
  dvector<int> nodes;
  int cs;
public:
  int size() { return go_template()->size_m; }
  int csum() { return go_template()->cs; }
  int dim() { return go_template()->dim_m; }
  int nperms() { return go_template()->perms_m; }
  const int* perm(int perm_indx) {
    return &(go_template()->perms_v.e(perm_indx,0)); 
  }
  GeomObject::Type type() { return go_template()->type; }
  GeomObjectSeq(const int *nodes_a) {
    int sz = size();
    nodes.set_chunk_size(sz);
    cs = 0;
    for (int j=0; j<sz; j++) {
      int node = nodes_a[j];
      cs += node;
      nodes.ref(j) = node;
    }
  }
};

class Edge : public GeomObjectSeq {
  Template *go_template() { return &GeomObjectSeq::EdgeTemplate; }
};
