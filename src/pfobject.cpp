//__INSERT_LICENSE__
// $Id: pfobject.cpp,v 1.5 2003/03/06 20:49:04 mstorti Exp $

#include <src/pfobject.h>
#include <src/texthash.h>
#include <petsc.h>

//#include <src/nullvort.h>

BasicObject_factory_t BasicObject_application_factory=NULL;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
BasicObject::~BasicObject() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
BasicObject *BasicObject::factory(string &type) {
  if (0) {} // tricky!!
  // Global (for all applications) BasicObject types
  // else if (type=="null_vort") return new null_vort; // now in NS
  // BasicObject types specific for this application
  else if (BasicObject_application_factory) {
    BasicObject *obj = (*BasicObject_application_factory)(type);
    if (obj) return obj;
  } else return NULL;
}

