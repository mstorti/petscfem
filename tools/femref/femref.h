// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.3 2004/11/17 23:45:41 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <src/generror.h>

typedef unsigned int Node;

class GeomObject {
public:
  virtual bool equal(GeomObject &)=0;
  virtual int dim()=0;
  static GeomObject *factory(int nnod,const int *nodes,
		      const char *type);
};

class OrientedEdge : public GeomObject {
private:
  int n1,n2;
  OrientedEdge(int n,const int *nodes) {
    if (n!=2) throw GenericError("OrientedEdge must have 2 nodes");
    n1 = nodes[0];
    n2 = nodes[1];
  };
public:
  int dim() { return 1; }
  bool equal(GeomObject &go) {
    OrientedEdge *edge2 = dynamic_cast<OrientedEdge *>(&go);
    if (!edge2) return false;
    return (n1==edge2->n1 && n2==edge2->n2);
  };
};

class Edge : public GeomObject {
private:
  int n1,n2;
  Edge(int n,const int *nodes) {
    if (n!=2) throw GenericError("Edge must have 2 nodes");
    n1 = nodes[0];
    n2 = nodes[1];
  };
public:
  int dim() { return 1; }
  bool equal(GeomObject &go) {
    Edge *edge2 = dynamic_cast<Edge *>(&go);
    if (!edge2) return false;
    return ((n1==edge2->n1 && n2==edge2->n2) ||
	    (n1==edge2->n2 && n2==edge2->n1));
  };
};

class OrientedTri : public GeomObject {
private:
  int nodes[3];
  OrientedTri(int n,const int *nodes_a) {
    if (n!=3) throw GenericError("OrientedTri must have 3 nodes");
    for (int j=0; j<3; j++) nodes[j] = nodes_a[j];
  };
public:
  int dim() { return 2; }
  bool equal(GeomObject &go) {
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

class Tri : public GeomObject {
private:
  int n1,n2,n3;
  int nodes[3];
  Tri(int n,const int *nodes_a) {
    if (n!=3) throw GenericError("Tri must have 3 nodes");
    for (int j=0; j<3; j++) nodes[j] = nodes_a[j];
  };
public:
  int dim() { return 2; }
  bool equal(GeomObject &go) {
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

#if 0
class ObjectTable

class NodeSet {
 public:
  /// Ctor 
  NodeSet(int ndim_a) : ndim(ndim_a), max(0) { }
  /** A node is identified by an index in `node_list'. 
      All indices between in range #[0,max-1)# that are not in 
      #invalid#, are nodes. **/
  class NodeHandle : public int { };
  /// Adds a new node with those specific coordinates
  NodeHandle add(double *coords);
  /// Remove that nodes
  void remove(NodeHandle node);
  
 private:
  /// The dimension of the space
  int ndim;
  /* Stores the coordinates.  A node is identified by an index in
     `node_list'.  All indices in range #[0,max-1)# that are not in
     #invalid#, are nodes. **/
  dvector<double> coords;
  /// Max node added so far
  int max;
  /// These are `holes' in the container
  set<int> invalid;
};
#endif

#if 0
class Mesh {
 public:
  class Object {};
  
 private:
  dvector<double> xnod();
  dvector<int> icone();
};
#endif

#endif
