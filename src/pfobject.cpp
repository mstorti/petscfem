//__INSERT_LICENSE__
// $Id: pfobject.cpp,v 1.3 2003/02/25 20:34:22 mstorti Exp $

#include <src/pfobject.h>
#include <src/texthash.h>
#include <petsc.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
BasicObject::~BasicObject() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
BasicObject *BasicObject::factory(string &type) {
  if (type=="null_vort") return new null_vort;
  else return NULL;
}

