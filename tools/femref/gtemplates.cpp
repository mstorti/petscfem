//__INSERT_LICENSE__
// $Id: gtemplates.cpp,v 1.9 2004/12/26 16:37:46 mstorti Exp $

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

int 
Tetra2TetraSplitter_rn[] 
= {0, 1,
   1, 2,
   0, 2,
   0, 3,
   1, 3,
   2, 3,GeomObject::NULL_NODE};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Tetra2TetraSplitterClass::
Tetra2TetraSplitterClass() {
  subobj_conn.mono(32);
  subobj_conn.reshape(2,8,4);
  ref_nodes.cat(Tetra2TetraSplitter_rn,
		GeomObject::NULL_NODE);
  ref_nodes.reshape(2,6,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Tetra2TetraSplitterClass::
size(GeomObject::Type t) const { 
  if (t==GeomObject::OrientedTetraT) return 8;
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
  t=GeomObject::OrientedTetraT;
  return &subobj_conn.e(j,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int 
Tetra2TetraSplitterClass::
nref_nodes() const { return 6; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
Tetra2TetraSplitterClass::
ref_node(int indx,
	 const GeomObject::Template *&tmpl,
	 int &nnod, 
	 const int *&nodes) {
  nnod = 2;
  nodes = &ref_nodes.e(indx,0);
  tmpl = &EdgeRefNodeTemplate;
}
