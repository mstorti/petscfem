// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.9 2004/11/19 01:33:41 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/linkgraph.h>
#include <src/generror.h>

typedef unsigned int Node;

class GeomObjectBasic {
public:
  virtual bool equal(GeomObjectBasic &)=0;
  virtual int dim()=0;
  static GeomObjectBasic *factory(int nnod,const int *nodes,
		      const char *type);
};

class GeomObject {
private:
  GeomObjectBasic *gobj;
public:
  GeomObject(GeomObjectBasic *go=NULL) : gobj(NULL) { }
  bool equal(GeomObject &go) { return gobj->equal(*go.gobj); }
  int dim() { return gobj->dim(); }
  void clear() { if (gobj) delete gobj; gobj=NULL; }
  ~GeomObject() { clear(); }
};

class OrientedEdge : public GeomObjectBasic {
private:
  int n1,n2;
  OrientedEdge(int n,const int *nodes) {
    if (n!=2) throw GenericError("OrientedEdge must have 2 nodes");
    n1 = nodes[0];
    n2 = nodes[1];
  };
public:
  int dim() { return 1; }
  bool equal(GeomObjectBasic &go) {
    OrientedEdge *edge2 = dynamic_cast<OrientedEdge *>(&go);
    if (!edge2) return false;
    return (n1==edge2->n1 && n2==edge2->n2);
  };
};

class Edge : public GeomObjectBasic {
private:
  int n1,n2;
  Edge(int n,const int *nodes) {
    if (n!=2) throw GenericError("Edge must have 2 nodes");
    n1 = nodes[0];
    n2 = nodes[1];
  };
public:
  int dim() { return 1; }
  bool equal(GeomObjectBasic &go) {
    Edge *edge2 = dynamic_cast<Edge *>(&go);
    if (!edge2) return false;
    return ((n1==edge2->n1 && n2==edge2->n2) ||
	    (n1==edge2->n2 && n2==edge2->n1));
  };
};

class OrientedTri : public GeomObjectBasic {
private:
  int nodes[3];
  OrientedTri(int n,const int *nodes_a) {
    if (n!=3) throw GenericError("OrientedTri must have 3 nodes");
    for (int j=0; j<3; j++) nodes[j] = nodes_a[j];
  };
public:
  int dim() { return 2; }
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

class Tri : public GeomObjectBasic {
private:
  int n1,n2,n3;
  int nodes[3];
  Tri(int n,const int *nodes_a) {
    if (n!=3) throw GenericError("Tri must have 3 nodes");
    for (int j=0; j<3; j++) nodes[j] = nodes_a[j];
  };
public:
  int dim() { return 2; }
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
  class GeomObjectAdaptor : public GeomObject{ 
    int dim() { }
    bool equal(GeomObjectBasic &go) { }
  };
  dvector<double> coords;
  dvector<int> tri;
  dvector<int> nod_adj;
  int ndim;
  LinkGraph lgraph;
  int nnod, nelem, nel;
public:
  class iterator { 
  private:
    GeomObjectAdaptor goa;
  public:
    GeomObject& operator*() { return goa; } 
    GeomObject* operator->() { return &goa; }
  };
  void get_adjacency(const GeomObject &g1,int dim,
		     std::list<iterator> &li) { 
    
  }
  Mesh(int ndim_a,int nel_a) 
    : ndim(ndim_a), nel(nel_a) { 
    coords.reshape(2,0,ndim);
    tri.reshape(2,0,nel);
  }
  void read(const char *node_file,const char *conn_file) {
    coords.cat(node_file);
    tri.cat(conn_file);
    nnod = coords.size(0);
    nel = tri.size(1);
    lgraph.init(nnod);
    nelem = tri.size(0);
    for (int ele=0; ele<nelem; ele++) {
      for (int k=0; k<nel; k++) {
	int node = tri.e(ele,k);
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
  }
};

#endif
