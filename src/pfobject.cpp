//__INSERT_LICENSE__
// $Id: pfobject.cpp,v 1.4 2003/02/26 00:54:50 mstorti Exp $

#include <src/pfobject.h>
#include <src/texthash.h>
#include <petsc.h>

#include <src/nullvort.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
BasicObject::~BasicObject() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
BasicObject *BasicObject::factory(string &type) {
  if (type=="null_vort") return new null_vort;
  else return NULL;
}

