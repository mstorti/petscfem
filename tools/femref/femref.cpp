//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.5 2004/11/20 21:16:27 mstorti Exp $

#include <string>

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

int GeomObject::Type GeomObject::type() {
  return gobj->type(); 
}

void GeomObject::clear() { 
  if (gobj) delete gobj; gobj=NULL; 
}

GeomObject::~GeomObject() { 
  clear(); 
}

int main() { 
  Mesh mesh(2,3);
  mesh.read("coord.dat","icone.dat");
}
