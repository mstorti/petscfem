//__INSERT_LICENSE__
//$Id: linkgraph.cpp,v 1.7 2002/07/24 01:19:54 mstorti Exp $

#include <src/linkgraph.h>

int LinkGraph::CHUNK_SIZE_DEF = 10000;
int LinkGraph::null = -1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraph::free_cell(int cell) {
  int_pair &fh = da.ref(M);
  int temp = fh.j;
  fh.j = cell;
  da.ref(cell).j = temp;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraph::erase(iterator q) {
  assert(q.graph == this);
  int p = da.ref(q.r).j, pp;
  while (p!=null) {
    pp = da.ref(p).j;
    free_cell(p);
    p = pp;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraph::erase(iterator first,iterator last) {
  for (iterator q=first; q!=last; q++) erase(q);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraph::init(int MM) {
  M=MM;
  da.resize(M+1);		// cell `M' is for the `free' cell list
  int_pair p(0,null);
  for (int j=0; j<=M; j++) da.ref(j) = p; // fill initial adjacency table
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int LinkGraph::available() {
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
void LinkGraph::list_insert(int header, int val) {
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
  Q.i = null;			// not needed?
  Q.j = null;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraph::set_ngbrs(int v,GSet &ngbrs) {
  int_pair *p = &da.ref(v);
  while (p->j != null) {
    ngbrs.insert(p->i);
    p = &da.ref(p->j);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int LinkGraph::size(int r) {
  int s=0;
  assert(s>=0 && s<M);
  int_pair *p = &da.ref(r);
  while (p->j != null) { s++; p = &da.ref(p->j); }
  return s;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Pack the row
void LinkGraph::pack(const Row & row,char *&buff) const {
  int n=row.size();
  BUFFER_PACK<int>(row.row,buff);
  BUFFER_PACK<int>(n,buff);
  Row::iterator q;
  for (q=row.begin(); q!=row.end(); q++)
    BUFFER_PACK<int>(*q,buff);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraph::unpack(Row & row,const char *&buff)  {
  int n,k;
  BUFFER_UNPACK<int>(row.row,buff);
  BUFFER_UNPACK<int>(n,buff);
  for (int j=0; j<n; j++) {
    BUFFER_UNPACK<int>(k,buff);
    row.insert(k);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// combine a row in the container
void LinkGraph::combine(const Row &row) {
  int j=row.row;
  Row::iterator q;
  for (q=row.begin(); q!=row.end(); q++) list_insert(j,*q);
}

