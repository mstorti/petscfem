//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.6 2004/11/21 20:21:53 mstorti Exp $

#include <string>
#include <limits.h>

using namespace std;

#include "./femref.h"

GeomObject::GeomObject() : gobj(NULL) { }

void GeomObject::init(int sz,const int *nodes,
		      GeomObject::Type t) { 
  clear();
#define CTOR(CLASS) \
else if(t==CLASS##T) gobj = new CLASS(sz,nodes)
  if (t==NULL_TYPE) assert(0);
  CTOR(OrientedEdge);
  CTOR(Edge);
  CTOR(OrientedTri);
  CTOR(Tri);
  else assert(0);
}

bool GeomObject::equal(GeomObject &go) { 
  return gobj->equal(*go.gobj); 
}

int GeomObject::dim() { 
  return gobj->dim(); 
}

int GeomObject::csum() { 
  return gobj->csum(); 
}

int GeomObject::size(GeomObject::Type t) {
  return gobj->size(t); 
}

int GeomObject::size() {
  return gobj->size(); 
}


int GeomObject::nperms() { 
  return gobj->nperms();
}

const int* GeomObject::perm(int perm_indx) {
  return gobj->perm(perm_indx);
}

GeomObject::Type GeomObject::type() {
  return gobj->type(); 
}

void GeomObject::clear() { 
  if (gobj) delete gobj; gobj=NULL; 
}
GeomObject::~GeomObject() { 
  clear(); 
}

GeomObject Mesh::id2obj(GeomObjectId id) {
  //  GeomObject 
}

int Edge::perms_v[2] = {1,0};

int GeomObjectSeq::NULL_NODE = INT_MAX;

GeomObjectSeq::GOTemplate
::GOTemplate(int sz,GeomObject::Type t,
	     int dim_a,int nperms_a,const int *perms) 
  : size_m(sz), type(t), dim_m(dim_a), nperms_m(nperms_a) {
  perms_v.resize(nperms_m*size_m);
  int *q = perms;
  for (int j=0; j<nperms_m*size_m; j++)
    perms_v.e(j) = *q++;
  perms_v.reshape(2,nperms_m,size_m);
}

static GeomObjectSeq::GOTemplate 
EdgeTemplate(2,GeomObject::EdgeT,1,1,{1,0}) { }

int main() { 
  Mesh mesh(2,3);
  mesh.read("coord.dat","icone.dat");
}

