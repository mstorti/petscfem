//__INSERT_LICENSE__
//$Id: linkgraph.cpp,v 1.3 2002/07/22 19:08:37 mstorti Exp $

#include <src/linkgraph.h>

int link_graph::CHUNK_SIZE_DEF = 10000;
int link_graph::null = -1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
link_graph::link_graph(int MM, const DofPartitioner *part,
		       MPI_Comm comm,int chunk_size=CHUNK_SIZE_DEF) 
  : M(MM), da(chunk_size) {
  if (M) init(MM);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void link_graph::init(int MM) {
  M=MM;
  da.resize(M+1);		// cell `M' is for the `free' cell list
  int_pair p(0,null);
  for (int j=0; j<=M; j++) da.ref(j) = p; // fill initial adjacency table
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int link_graph::available() {
  int_pair &fh = da.ref(M);
  if (fh.j == null) {
    da.push(int_pair());
    return da.size()-1;
  } else {
    int free = fh.j;
    fh.j = da.ref(fh.j).j;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void link_graph::list_insert(int header, int val) {
  int_pair *p = &da.ref(header);
  while (p->j != null) {
    if (p->i == val) return;
    p = &da.ref(p->j);
  }
  int q = available();
  assert(q!=null);
  p->j = q;
  p->i = val;
  int_pair &Q = da.ref(q);
  Q.i = null;			// not needed
  Q.j = null;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void link_graph::set_ngbrs(int v,GSet &ngbrs) {
  int_pair *p = &da.ref(v);
  while (p->j != null) {
    ngbrs.insert(p->i);
    p = &da.ref(p->j);
  }
}
