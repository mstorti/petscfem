// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.25 2004/11/23 13:53:58 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/generror.h>
#include <src/linkgraph.h>

/** Represents a geomtric object `outside' its container.
    Basically, it represents a sequence of nodes and a 
    topological shape (template). */ 
class GeomObject {
public:
  /** Types of objects may be identified by a #Type#
      value, or as a pointer to a #Template# class. 
   */ 
  enum Type { NULL_TYPE=0, 
	      OrientedEdgeT, EdgeT, 
	      OrientedTriT, TriT,
	      OrientedTetraT, TetraT};
  /// This is a value that a real node can't take. 
  static int NULL_NODE;
  /// Contains information about a shape
  class Template {
  public:
    // friend class GeomObject;
    // friend class Mesh;
    /// Contains the permutations that leave the shape invariant
    dvector<int> perms_v;
    /// Number of nodes, dimension, number of permutations
    int size_m,dim_m,nperms_m;
    /// The corresponding enum #Type#. 
    GeomObject::Type type;
    /// The name of the shape
    const char *label;
    /// Ctor from info
    Template(int sz,GeomObject::Type t,
	     int dim_a,int nperms_a,
	     const int *perms,const char *label_a);
    /// The number of objects of type #t# this shape has. 
    virtual int size(Type t) const { return 0; }
    /// Local nodes connected to subobject #j# of type #t#. 
    virtual const int* nodes(Type t,int j) const { assert(0); }
  };
  /// Ptr to shape structure
  Template *go_template;
  /** This function returns the shape ptr for a given
      enum #Type# */ 
  static Template* get_template(Type t);
private:
  /// The check-sum for this object, has been converted to canonical?
  int cs, canonical;
  /// The nodes for this object
  dvector<int> nodes_m;
public:
  /// Number of nodes for this object
  int size() const { return go_template->size_m; }
  /// The number of objects of type #t# this shape has. 
  int size(Type t) const { return go_template->size(t); }
  /// The checksum for this object
  int csum() const { return cs; }
  /// Dimension of the object
  int dim() const { return go_template->dim_m; }
  /// Number of permutations that leave this object invariant. 
  int nperms() const { return go_template->nperms_m; }
  /// Reset to the empty object
  void clear() { go_template=NULL; nodes_m.clear(); canonical=0; cs=0; }
  /// Set to an object with the specific info
  void init(Type t,const int *nodes_a=NULL);
  /// Local nodes for #j# permutation 
  const int* perm(int j) const {
    return &(go_template->perms_v.e(j,0)); 
  }
  /// Returns the #enum# type. 
  Type type() const { return go_template->type;}
  /// Constructor of empty object
  GeomObject() : canonical(0), go_template(NULL) { }
  /// Constructor with info
  GeomObject(Type t,const int *nodes_a=NULL);
  /** Bring to canonical ordering by permuting the
      nodes.  The canonical ordering is that
      permutation of the nodes that gives the lowest
      lexicographical order of the sequence of nodes.
      The sequence of nodes in canonical order is
      a unique identification of the object. 
   */ 
  void make_canonical();
  /// Compare two objects
  bool equal(GeomObject &go);
  /// Finds the subobject that matches object #go#. 
  int find(GeomObject &go) const;
  /// Prints the object
  void print(const char*s = NULL) const;
  /// The nodes of this object
  const int* nodes() const { return nodes_m.buff(); }
  /** Set #go# to be the #j#-th subobject of type #t# in
      #*this# */ 
  void set(Type t,int j,GeomObject &go) const;
};

/** A mesh is a container of geometrical objects
    linked (as in a graph).  Two objects are linked
    if they share a node. For reasons of efficiency,
    meshes are not stored by storing the whole graph. */ 
class Mesh {
public:
  /** An iterator is a position in the
      graph. Basically it is the number of element
      in the mesh plus the number that identifies
      the subobject inside the element. Iterators do
      not represent uniquely the object,
      i.e. several iterators can point to the same
      object. */ 
  class iterator {
  public:
    // friend class Mesh;
    /// The big object number in the mesh
    int obj;
    /** The type of the geometric object
	identified by this iterator */ 
    GeomObject::Type t;
    /// The position of the subobject in the large object
    int subobj;
    /// Ctor from the data
    iterator(int obj_a,GeomObject::Type ta,int subobj_a)
      : obj(obj_a), t(ta), subobj(subobj_a) { }
    /// Ctor of empty iterator
    iterator() : obj(-1), t(GeomObject::NULL_TYPE), 
		 subobj(-1) { }
    /// Set to data 
    void set(int obj_a,GeomObject::Type ta,int subobj_a) {
      obj = obj_a; t = ta; subobj = subobj_a; 
    }
  };
  /** Set object #go# to the object pointed
      by iterator #it# */
  virtual void set(iterator it,GeomObject &go)=0;
  /** Returns one iterator that points to the 
      specified object. Remember that may be many. */ 
  virtual iterator find(GeomObject &go)=0;
  /** Finds the list of all iterators that point to
      object #go#. */ 
  virtual void find(GeomObject &go,list<iterator> &its)=0;
  /** Get list #adj# of iterators of shape #t# that
      that are connected to iterator #t#  */ 
  virtual void get_adjacency(iterator it,GeomObject::Type t,
			     list<iterator> &adj)=0;
  /// Determine whether this iterator is empty. 
  virtual bool is_end(iterator it)=0;
};

/** This is a container of objects formed by 
    all large objects of shape #tmpl#. */ 
class UniformMesh : public Mesh {
private:
  /// The coordinates of the nodes
  dvector<double> coords;
  /** The connectivities of the mesh. List of
      nodes of each large object. */
  dvector<int> connec;
  /** Adjacency list node-to-element. Stores the list of
      elements connected to a given node. */ 
  dvector<int> n2e;
  /** Pointers to ranges in #n2e#, elements connected 
      to node #j# are in range #[n2e_ptr[j],n2e_ptr[j+1])#. 
      (#n2e_ptr# is size #nnod+1#). */
  dvector<int> n2e_ptr;
  /// The dimension of the space
  int ndim;
  /** Number of nodes, number of elements, number of
      nodes connected to an element */
  int nnod, nelem, nel;
  
  const GeomObject::Template *tmpl;
  void find(GeomObject &go,list<iterator> *its,iterator &it);

public:
  UniformMesh(GeomObject::Template &tmpl_a,int ndim_a) 
    : tmpl(&tmpl_a), nel(tmpl_a.size_m), ndim(ndim_a) { 
    coords.reshape(2,0,ndim);
    connec.reshape(2,0,nel);
  }

  void read(const char *node_file,const char *conn_file) {
    LinkGraph lgraph;
    coords.cat(node_file).defrag();
    connec.cat(conn_file).defrag();
    nnod = coords.size(0);
    nel = connec.size(1);
    lgraph.init(nnod);
    nelem = connec.size(0);
    for (int ele=0; ele<nelem; ele++) {
      for (int k=0; k<nel; k++) {
	int node = connec.e(ele,k);
	printf("add node %d, elem %d\n",node,ele);
	lgraph.add(node,ele);
      }
    }
    GSet ngbrs;
    for (int node=0; node<nnod; node++) {
      printf("node %d, elems ",node);
      ngbrs.clear();
      lgraph.set_ngbrs(node,ngbrs);
      GSet::iterator p = ngbrs.begin();
      while (p!=ngbrs.end()) printf("%d ",*p++);
      printf("\n");
    }

    // Pass the connectivity in graph `lgraph' to
    // node_list
    n2e_ptr.resize(nnod+1);
    for (int j=0; j<nnod; j++) n2e_ptr.ref(j) = 0;
    int nadj = 0;
    for (int nod=0; nod<nnod; nod++) {
      n2e_ptr.e(nod) = nadj;
      ngbrs.clear();
      lgraph.set_ngbrs(nod,ngbrs);
      nadj += ngbrs.size();
      GSet::iterator p = ngbrs.begin();
      while (p!=ngbrs.end()) n2e.push(*p++);
    }
    n2e_ptr.ref(nnod) = nadj;
    lgraph.clear();
  }
  
  void set(iterator it,GeomObject &go);
  iterator find(GeomObject &go);
  void find(GeomObject &go,list<iterator> &its);
  void get_adjacency(iterator it,GeomObject::Type t,
		     list<iterator> &adj) { }
  bool is_end(iterator it) { 
    return it.obj<0 || it.t==GeomObject::NULL_TYPE 
      || it.subobj<0; 
  }
};

#endif
