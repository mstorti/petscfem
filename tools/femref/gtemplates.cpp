//__INSERT_LICENSE__
// $Id: gtemplates.cpp,v 1.1 2004/12/05 20:04:47 mstorti Exp $

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
  subobj_conn.a_resize(2,8,4);
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
