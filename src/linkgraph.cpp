//__INSERT_LICENSE__
//$Id: linkgraph.cpp,v 1.2 2002/07/22 15:45:12 mstorti Exp $

#include <src/linkgraph.h>

int link_graph::CHUNK_SIZE_DEF = 10000;
int link_graph::null = -1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
link_graph::link_graph(int MM, const DofPartitioner *part,
		       MPI_Comm comm,int chunk_size=CHUNK_SIZE_DEF) 
  : M(MM), da(chunk_size) {
  da.resize(M+1); //cell `M' is for the `free' cell list
  int_pair p(0,null);
  for (int j=0; j<=M; j++) da.ref(j) = p;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int link_graph::available() {
  int_pair &fh = da.ref(M);
  if (fh.j==-1) {
    da.push(int_pair());
    return size();
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
  da.ref(q).i = val;
}
