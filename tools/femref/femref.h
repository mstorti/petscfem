// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.12 2004/11/21 14:51:54 mstorti Exp $
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
  static int size(Type t);
  /// Number of nodes adjacent to this object. 
  static int size();
  /// Type of this object. 
  int GeomObject::Type type();
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
  virtual static int size();
  /** Returns the number of object of type #t#
      connected to #*this#. */ 
  virtual static 
  int size(GeomObject::Type t) {return 0; }
  
  virtual GeomObject::Type type() { 
    return GeomObject::NULL_TYPE; 
  }
  virtual int nperms() { return 0; }
  virtual const int perm(int perm_indx);
  virtual const int *
  static GeomObjectBasic 
  *factory(int sz,const int *nodes, GeomObject::Type t);
};

class OrientedEdge : public GeomObjectBasic {
private:
  int n1,n2,cs;
public:
  friend class Edge;
  OrientedEdge(int n,const int *nodes) {
    if (n!=2) throw GenericError("OrientedEdge must have 2 nodes");
    n1 = nodes[0];
    n2 = nodes[1];
    cs = n1 + n2;
  };
  int dim() { return 1; }
  int csum() { return cs; }
  static int size(GeomObject::Type t) { assert(0); }
  static int size() { return 2; }
  GeomObject::Type type() { 
    return GeomObject::OrientedEdgeT; 
  }
  bool equal(GeomObjectBasic &go) {
    OrientedEdge *edge2 = dynamic_cast<OrientedEdge *>(&go);
    if (!edge2) return false;
    return (n1==edge2->n1 && n2==edge2->n2);
  };
};

class Edge : public OrientedEdge {
public:
  Edge(int n,const int *nodes) : OrientedEdge(n,nodes) { }
  bool equal(GeomObjectBasic &go) {
    Edge *edge2 = dynamic_cast<Edge *>(&go);
    if (!edge2) return false;
    return ((n1==edge2->n1 && n2==edge2->n2) ||
	    (n1==edge2->n2 && n2==edge2->n1));
  };
  GeomObject::Type type() { 
    return GeomObject::EdgeT; 
  }
};

class OrientedTri : public GeomObjectBasic {
private:
  int nodes[3],cs;
public:
  friend class Tri;
  OrientedTri(int n,const int *nodes_a) {
    if (n!=3) throw GenericError("OrientedTri must have 3 nodes");
    cs = 0;
    for (int j=0; j<3; j++) {
      nodes[j] = nodes_a[j];
      cs += nodes[j];
    }
  };
  int dim() { return 2; }
  int csum() { return cs; }
  static int size() { return 3; }
  static int size(GeomObject::Type t) { 
    if (t==EdgeT) return 3; 
    else assert(0);
  }
  GeomObject::Type type() { 
    return GeomObject::OrientedTriT; 
  }
  bool equal(GeomObjectBasic &go) {
    OrientedTri *tri2 = dynamic_cast<OrientedTri *>(&go);
    if (!tri2) return false;
    int n0 = tri2->nodes[0];
    int j;
    for (j=0; j<3; j++) if (nodes[j] = n0) break;
    if (j==3) return false;
    for (int k=1; k<3; k++)
      if (tri2->nodes[k] != nodes[(j+k)%3]) return false;
    return true;
  };
};

class Tri : public OrientedTri {
public:
  Tri(int n,const int *nodes) : OrientedTri(n,nodes) { }
  static int size(GeomObject::Type t) { assert(0); }
  GeomObject::Type type() { 
    return GeomObject::TriT; 
  }
  bool equal(GeomObjectBasic &go) {
    Tri *tri2 = dynamic_cast<Tri *>(&go);
    if (!tri2) return false;
    int n0 = tri2->nodes[0];
    int j;
    for (j=0; j<3; j++) 
      if (nodes[j] = n0) break;
    if (j==3) return false;
    bool flag = true;
    // Check for same orientation
    for (int k=1; k<3; k++) {
      if (tri2->nodes[k] != nodes[(j+k)%3]) {
	flag = false; break;
      }
    }
    if (flag) return true;
    // Check for reverse orientation
    for (int k=1; k<3; k++) {
      if (tri2->nodes[k] != nodes[(j-k)%3]) {
	flag = false; break;
      }
    }
    return flag;
  };
};

class Mesh {
private:
  class GeomObjectId {
  private:
    int elem, local_number;
    GeomObject::Type t;
    GeomObjectId(int elem_a,int lna,
		 GeomObject::Type ta) 
      : elem(elem_a), local_number(lna), t(ta) {};
  public:
    friend class Mesh;
    int csum() { };
    GeomObject operator*();
  };
  dvector<double> coords;
  dvector<int> connec;
  dvector<int> n2e;
  dvector<int> n2e_ptr;
  int ndim;
  LinkGraph lgraph;
  int nnod, nelem, nel;
public:
#if 0
  class iterator { 
  private:
    GeomObjectAdaptor goa;
  public:
    GeomObject& operator*() { return goa; } 
    GeomObject* operator->() { return &goa; }
  };
#endif
  GeomObject id2obj(GeomObjectId id);
  void get_adjacency(const GeomObject &g1,int dim,
		     std::list<GeomObjectId> &li) { }
  Mesh(int ndim_a,int nel_a) 
    : ndim(ndim_a), nel(nel_a) { 
    coords.reshape(2,0,ndim);
    connec..reshape(2,0,nel);
  }
  void read(const char *node_file,const char *conn_file) {
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
    for (int j=0; j<nnod; j++) n2e.ref(j) = 0;
    int nadj = 0;
    for (int nod=0; nod<nnod; nod++) {
      n2e_ptr.e(nod) = nadj;
      ngbrs.clear();
      lgraph.set_ngbrs(nod,ngbrs);
      nadj += ngbrs.size();
      GSet::iterator p = ngbrs.begin();
      while (p!=ngbrs.end()) n2e.push(*p);
    }
    n2e_ptr.ref(nnod) = nadj;
    lgraph.clear();
  }

  void list_faces() {
    GeomObject elem;
    for (int e=0; e<nelem; e++) {
      elem.init(e,connec.e(e,0),OrientedTriT);
      int nedge = elem.size(Edge);
      for (int le=0; le<nedge; le++) {
	GeomObjectId lid(e,le,EdgeT);
	int nodes[2];
	nodes[0] = connec.e(ele,0);
	nodes[1] = connec.e(ele,1);
	Edge edge(2,nodes);
	int csum = nodes[1] + nodes[1];
	int p0 = n2e_ptr.ref(nodes[0]);
	int p0e = n2e_ptr.ref(nodes[0]+1);
	int p1 = n2e_ptr.ref(nodes[1]);
	int p1e = n2e_ptr.ref(nodes[1]+1);
	while (p0<p0e && p1<p1e) {
	  int ele0 = n2e.ref(p0);
	  int ele1 = n2e.ref(p1);
	  if (ele0 && ele1) {
	  } else if (ele0<ele1) ele0++;
	  else ele1++;
	}
      }
    }
  }
};

#endif
