//__INSERT_LICENSE__
// $Id: gtemplates.cpp,v 1.2 2004/12/05 21:43:48 mstorti Exp $

#include <string>
#include <list>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int 
Tetra2TetraSplitter_v[] 
= {4,7,6,0,
   4,5,8,1,
   5,6,9,2,
   8,9,7,3,
   4,6,7,8,
   4,5,6,8,
   5,9,6,8,
   9,7,8,6,GeomObject::NULL_NODE};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Tetra2TetraSplitterClass::
Tetra2TetraSplitterClass() {
  subobj_conn.a_resize(2,8,4,2);
  int nso = 8;
  int nel = 4;
  int nelso = 4;
  for (int j=0; j<nso; j++) {
  }
  subobj_conn.set(Tetra2TetraSplitter_v);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Tetra2TetraSplitterClass::
size(GeomObject::Type t) const { 
  if (t==GeomObject::TetraT) return 8;
  else return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const int* 
Tetra2TetraSplitterClass::
nodes(GeomObject::Type t,int j) { 
  if (t==GeomObject::TetraT) 
    return &subobj_conn.e(j,0);
  else assert(0);
}
