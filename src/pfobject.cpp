//__INSERT_LICENSE__
// $Id: pfobject.cpp,v 1.10 2005/04/09 21:45:35 mstorti Exp $

#include <src/pfobject.h>
#include <src/texthash.h>
#include <petsc.h>

#include <src/interpola.h>

BasicObject_factory_t *BasicObject_application_factory=NULL;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
BasicObject::~BasicObject() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
BasicObject *BasicObject::factory(string &type) {
  if (0) {} // tricky!!
  // else if (type=="interpolation") return new interpolation;

  // BasicObject types specific for this application
  if (BasicObject_application_factory) {
    BasicObject *obj = (*BasicObject_application_factory)(type);
    if (obj) return obj;
  } 

  return NULL;
}

