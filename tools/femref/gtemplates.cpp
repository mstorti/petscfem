//__INSERT_LICENSE__
// $Id: gtemplates.cpp,v 1.6 2004/12/14 18:05:54 mstorti Exp $

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
  int nel = 4;
  int nso = 8;
  int nelso = 4;
  for (int j=0; j<nso; j++) {
    for (int k=0; k<nelso; k++) {
      int n2=0, n1 = Tetra2TetraSplitter_v[nelso*j+k];
      if (n1>=nel) switch(n1) {
      case 4: n1=0; n2=1; break;
      case 5: n1=1; n2=2; break;
      case 6: n1=0; n2=2; break;
      case 7: n1=0; n2=3; break;
      case 8: n1=1; n2=3; break;
      case 9: n1=2; n2=3; break;
      default: assert(0);
      }
      subobj_conn.e(j,k,0) = n1;
      subobj_conn.e(j,k,1) = n2;
    }
  }
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
nodes(GeomObject::Type t,int j) const { 
  if (t==GeomObject::TetraT) 
    return &subobj_conn.e(j,0);
  else assert(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Tetra2TetraSplitterClass::
size() const { return 8; };

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const int* 
Tetra2TetraSplitterClass::
nodes(int j,GeomObject::Type &t) const { 
  t=GeomObject::TetraT;
  return &subobj_conn.e(j,0);
}
