//__INSERT_LICENSE__
// $Id: pfobject.cpp,v 1.2 2003/02/25 13:31:41 mstorti Exp $

#include <src/pfobject.h>
#include <src/texthash.h>
#include <petsc.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
BasicObject::~BasicObject() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class dummy_obj : public BasicObject {
private:
  TextHashTable thash;
public:
  ~dummy_obj() {}
  void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) { 
    thash.read(fstack);
    thash.print("Read dummy obj:");
    PetscFinalize();
    exit(0);
  }
};
	
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
BasicObject *BasicObject::factory(string &type) {
  if (type=="dummy_obj") return new dummy_obj;
  else return NULL;
}

