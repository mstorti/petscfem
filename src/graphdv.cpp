//__INSERT_LICENSE__
//$Id: graphdv.cpp,v 1.2 2002/07/22 15:45:12 mstorti Exp $

#include <src/graphdv.h>

int graphdv::CHUNK_SIZE_DEF = 10000;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graphdv::resync() { 
  if (modif) {
    // printf("resyncing at size: %d\n",da.size());
    da.sort();
    da.shrink(da.remove_unique());
    modif = 0;
    // printf("end resync.\n",da.size());
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
graphdv::graphdv(int chunk_size=CHUNK_SIZE_DEF) : da (chunk_size) { 
  ordered = 0; MAX_INIT = 1000; 
  max=MAX_INIT; modif=0; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graphdv::add(int j,int k) {
  // construct pair
  int_pair p(j,k);
  // insert if not in ordered region
  if (!da.find(p,0,ordered)) {
    // insert at the back
    da.push(p); modif=1; 
    // if size greater than threshold resync
    int ds = da.size();
    if (ds > max) {
      // launch automatic resync
      resync();
      int new_max = 2 * da.size();
      if (new_max > max) {
	max = new_max;
      }
    }
  } 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graphdv::print(const char *s=NULL) {
  resync(); 
  if (s) printf("%s",s);
  int ss = da.size();
  for (int j=0; j<ss; j++) {
    int_pair &q = da.ref(j);
    printf("(%d,%d) ",q.i,q.j);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graphdv::set_ngbrs(int v,vector<int> &ngbrs) {
  resync();
  int f,e,q;
  f = da.bsearch(int_pair(v,0));
  e = da.bsearch(int_pair(v+1,0));
  for (q=f; q<e; q++) ngbrs.push_back(da.ref(q).j);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graphdv::set_ngbrs(int v,GSet &ngbrs) {
  resync();
  int f,e,q;
  f = da.bsearch(int_pair(v,0));
  e = da.bsearch(int_pair(v+1,0));
  for (q=f; q<e; q++) ngbrs.insert(da.ref(q).j);
}
