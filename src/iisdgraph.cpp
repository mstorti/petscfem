//__INSERT_LICENSE__
//$Id: iisdgraph.cpp,v 1.1.2.5 2001/12/18 12:59:49 mstorti Exp $

// fixme:= this may not work in all applications
int MY_RANK,SIZE;

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
void StoreGraph::set_ngbrs(int j,set<int> &ngbrs_v) {
  GMap::iterator q;
  q = lgraph.find(j);
  if (q != lgraph.end()) ngbrs_v = q->second;
}

int DGMap::size_of_pack(const GRow &q) const {
  /// Contains #int k# and the set of neighobrs #set<int> s#. We store
  //#n+2# where #n=s.size()# integers in the sequence: #k, n, s[0],
  //... s[n-1]#.
  return sizeof(int)*(2+q.second.size());
}

void DGMap::pack(const GRow &p,char *&buff) const {
  int n;
  set<int>::const_iterator q,qe;
  const set<int> &w = p.second;
  BUFFER_PACK(p.first,buff);
  n = w.size();
  BUFFER_PACK(n,buff);
  qe=w.end();
  for (q=w.begin(); q!=qe; q++) {
    BUFFER_PACK(*q,buff);
#ifdef DEBUG_SG  // debug:=
    printf("[%d] packing %d -> %d\n",MY_RANK,p.first,*q);
#endif
  }  
}

void DGMap::unpack(GRow &p,const char *&buff) {
  set<int> &w = p.second;
  int k,n,q;

  BUFFER_UNPACK(p.first,buff);
  BUFFER_UNPACK(n,buff);
  
  for (k=0; k<n; k++) {
    BUFFER_UNPACK(q,buff);
#ifdef DEBUG_SG  // debug:=
    printf("[%d] unpacking %d -> %d\n",MY_RANK,p.first,q);
#endif
    w.insert(q);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DGMap::combine(const GRow &p) {
  GMap::iterator iter;
  set<int>::const_iterator q,qe;
  iter = find(p.first);
  if (iter != end()) {
    const set<int> &r = p.second;
    qe = r.end();
    for (q=r.begin(); q!=qe; q++) {
#ifdef DEBUG_SG  // debug:=
      printf("[%d] adding %d -> %d\n",MY_RANK,p.first,*q);
#endif
      iter->second.insert(*q);
    }
  } else {
#ifdef DEBUG_SG // debug:=
    const set<int> &r = p.second;
    qe = r.end();
    for (q=r.begin(); q!=qe; q++) {
      printf("[%d] adding %d -> %d\n",MY_RANK,p.first,*q);
    }
#endif
    insert(p);
  }
};

void StoreGraph::print() const {
  int p;
  DGMap::const_iterator q,qe;
  set<int>::const_iterator s,se;

  qe = lgraph.end();
  for (q = lgraph.begin(); q!=qe; q++) {
    const set<int> &w = q->second;
    PetscSynchronizedPrintf(comm,"[%d] %d -> ",MY_RANK,q->first);
    se = w.end();
    for (s=w.begin(); s!=se; s++)
      PetscSynchronizedPrintf(comm,"%d ",*s);
    PetscSynchronizedPrintf(comm,"\n");
  }
  PetscSynchronizedFlush(comm);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
#if 0
void StoreGraph::set_ngbrs(int loc1,set<int> &ngbrs_v) {
  int pos,loc2,dof2;
  pos = loc2dof[loc1]+k1;
  while (1) {
    nodep = (Node *)da_ref(da,pos);
    if (nodep->next==-1) break;
    // loc2:= number of global dof connected to `'
    dof2 = nodep->val;
    if (k1<=dof2 && dof2<=k2 && !flag[dof2] ) {
      loc2 = dof2loc[dof2-k1];
      ngbrs_v.insert(loc2);
    }
    pos = nodep->next;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double StoreGraph::weight(int elem) {
  return 1.;
}
#endif
