// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.20 2004/11/22 23:13:06 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/linkgraph.h>
#include <src/generror.h>

class GeomObject {
public:
  /// Types of objects
  enum Type { NULL_TYPE=0, 
	      OrientedEdgeT, EdgeT, 
	      OrientedTriT, TriT,
	      OrientedTetraT, TetraT};
  static int NULL_NODE;
  class Template {
  public:
    friend class GeomObject;
    friend class Mesh;
    dvector<int> perms_v;
    int size_m,dim_m,nperms_m;
    GeomObject::Type type;
    const char *label;
    Template(int sz,GeomObject::Type t,
	     int dim_a,int nperms_a,
	     const int *perms,const char *label_a);
    virtual int size(Type t) const { return 0; }
    virtual const int* nodes(Type t,int j) const { assert(0); }
  };
  Template *go_template;
  static Template* get_template(Type t);
private:
  int cs, canonical;
  dvector<int> nodes_m;
public:
  int size() const { return go_template->size_m; }
  int size(Type t) const;
  int csum() const { return cs; }
  int dim() const { return go_template->dim_m; }
  int nperms() const { return go_template->nperms_m; }
  void clear() { go_template=NULL; nodes_m.clear(); canonical=0; cs=0; }
  void init(Type t,const int *nodes_a=NULL);
  const int* perm(int perm_indx) const {
    return &(go_template->perms_v.e(perm_indx,0)); 
  }
  Type type() const { return go_template->type;}
  GeomObject() : canonical(0), go_template(NULL) { }
  GeomObject(Type t,const int *nodes_a=NULL);
  void make_canonical();
  /// Compare two objects
  bool equal(GeomObject &go);
  void print(const char*s = NULL) const;
  const int* nodes() const { return nodes_m.buff(); }
};

class Mesh {
public:
  class iterator {
  public:
    friend class Mesh;
    int obj;
    GeomObject::Type t;
    int subobj;
    iterator(int obj_a,GeomObject::Type ta,int subobj_a)
      : obj(obj_a), t(ta), subobj(subobj_a) { }
  };
  virtual void set(iterator it,GeomObject &go)=0;
  virtual iterator find(const GeomObject &go)=0;
  virtual void get_adjacency(iterator it,GeomObject::Type t,
			     list<iterator> &adj)=0;
};

class UniformMesh : public Mesh {
private:
  dvector<double> coords;
  dvector<int> connec;
  dvector<int> n2e;
  dvector<int> n2e_ptr;
  int ndim;
  LinkGraph lgraph;
  int nnod, nelem, nel;
  const GeomObject::Template *tmpl;

public:
  UniformMesh(GeomObject::Template &tmpl_a,int ndim_a) 
    : tmpl(&tmpl_a), nel(tmpl_a.size_m), ndim(ndim_a) { 
    coords.reshape(2,0,ndim);
    connec.reshape(2,0,nel);
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
  iterator find(const GeomObject &go);
  void get_adjacency(iterator it,GeomObject::Type t,
		     list<iterator> &adj) { }
};

#endif
