//__INSERT_LICENSE__
//$Id: iisdgraph.cpp,v 1.1.2.3 2001/12/17 00:03:57 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

#include <map>
#include <typeinfo>

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
  q=w.end();
  for (q=w.begin(); q!=qe; q++) 
    BUFFER_PACK(*q,buff);
}

void DGMap::unpack(GRow &p,const char *&buff) {
  set<int> &w = p.second;
  int k,n,q;

  BUFFER_UNPACK(p.first,buff);
  BUFFER_UNPACK(n,buff);
  
  for (k=0; k<n; k++) {
    BUFFER_UNPACK(q,buff);
    w.insert(q);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DGMap::combine(const GRow &p) {
  GMap::iterator iter;
  set<int>::const_iterator q,qe;
  iter = find(p.first);
  if (iter != end()) {
    const set<int> &insrow = p.second;
    qe = insrow.end();
    for (q=insrow.begin(); q!=qe; q++) iter->second.insert(*q);
  } else {
    insert(p);
  }
};

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
