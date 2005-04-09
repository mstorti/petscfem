//__INSERT_LICENSE__
// $Id: pfobject.cpp,v 1.9 2005/04/09 21:38:26 mstorti Exp $

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
  else if (BasicObject_application_factory) {
    return (*BasicObject_application_factory)(type);
  } else return NULL;
}

