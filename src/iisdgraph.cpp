//__INSERT_LICENSE__
//$Id: iisdgraph.cpp,v 1.5 2002/07/21 22:38:30 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

//#define DEBUG_SG  // debug:=
//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

#include <map>
#include <typeinfo>

#include <petsc.h>

#include <src/graph.h>
#include <src/distcont.h>
#include <src/iisdgraph.h>
#include <src/buffpack.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void StoreGraph1::set_ngbrs(int j,GSet &ngbrs_v) {
  GMap::iterator q;
  q = lgraph.find(j);
  if (q != lgraph.end()) ngbrs_v = q->second;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DGMap::size_of_pack(const GRow &q) const {
  /// Contains #int k# and the set of neighobrs #set<int> s#. We store
  //#n+2# where #n=s.size()# integers in the sequence: #k, n, s[0],
  //... s[n-1]#.
  return sizeof(int)*(2+q.second.size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DGMap::pack(const GRow &p,char *&buff) const {
  // number of rows in the map
  int n;
  GSet::const_iterator q,qe;
  // the set of integers pointed to p.first
  const GSet &w = p.second;
  // Pack the row number
  BUFFER_PACK(p.first,buff);
  n = w.size();
  // pack the number of items in the set of neighbors
  BUFFER_PACK(n,buff);
  // loop over the neighbors
  qe=w.end();
  for (q=w.begin(); q!=qe; q++) {
    // pack the neighbor
    BUFFER_PACK(*q,buff);
#ifdef DEBUG_SG  // debug:=
    printf("[%d] packing %d -> %d\n",MY_RANK,p.first,*q);
#endif
  }  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DGMap::unpack(GRow &p,const char *&buff) {
  // The et of neighbors
  GSet &w = p.second;
  int k,n,q;

  // unpack the row number
  BUFFER_UNPACK(p.first,buff);
  // unpack the number of ngbrs in the set
  BUFFER_UNPACK(n,buff);
  
  // loop over the ngbrs
  for (k=0; k<n; k++) {
    // unpack a ngbr
    BUFFER_UNPACK(q,buff);
#ifdef DEBUG_SG  // debug:=
    printf("[%d] unpacking %d -> %d\n",MY_RANK,p.first,q);
#endif
    // insert it in the merged set
    w.insert(q);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DGMap::combine(const GRow &p) {
  GMap::iterator iter;
  GSet::const_iterator q,qe;
  // Look if the row number is already in the graph
  iter = find(p.first);
  if (iter != end()) {
    // set of integers for this row
    const GSet &r = p.second;
    // loop over the set of ngbr to merge
    qe = r.end();
    for (q=r.begin(); q!=qe; q++) {
#ifdef DEBUG_SG  // debug:=
      printf("[%d] adding %d -> %d\n",MY_RANK,p.first,*q);
#endif
      // insert new neighbor
      iter->second.insert(*q);
    }
  } else {
#ifdef DEBUG_SG // debug:=
    const GSet &r = p.second;
    qe = r.end();
    for (q=r.begin(); q!=qe; q++) {
      printf("[%d] adding %d -> %d\n",MY_RANK,p.first,*q);
    }
#endif
    // insert the whole pair <row, set of ngbrs>
    insert(p);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void StoreGraph1::print() const {
  // print the map using `PetscSynchronizedPrintf'
  int p;
  DGMap::const_iterator q,qe;
  GSet::const_iterator s,se;

  // loop over rows
  qe = lgraph.end();
  for (q = lgraph.begin(); q!=qe; q++) {
    const GSet &w = q->second;
    // print processor and row number
    PetscSynchronizedPrintf(comm,"[%d] %d -> ",MY_RANK,q->first);
    // loop over ngbrs
    se = w.end();
    for (s=w.begin(); s!=se; s++)
      // print ngbr
      PetscSynchronizedPrintf(comm,"%d ",*s);
    PetscSynchronizedPrintf(comm,"\n");
  }
  PetscSynchronizedFlush(comm);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void StoreGraph1::clear() {
  Graph::clear();
  lgraph.clear();
}
