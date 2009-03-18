///__INSERT_LICENSE__
//$Id mstorti-v6-branch-1.0.1-5-ge86f38c Wed Sep 19 13:06:15 2007 -0300$

#include <cmath>
#include <cstdio>

using namespace std;
#include "fem.h"
#include "fastmat2.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2
::CacheCtx2::CacheCtx2() :
  branch_indx(-1), 
  branch_p(NULL) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx2
::Branch::Branch() : indx(-1) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2::jump(Branch &b) {
  if (b.indx==-1) {
    branchv.push_back(clist_t());
    b.indx = branchv.size()-1;
  }
  branch_indx = b.indx;
  branch_p = &branchv[branch_indx];
  q = branch_p->begin();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::clear() { branchv.clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCache * 
FastMat2::CacheCtx2::step() {
  FastMatCache *cache=NULL;
  if (q != branch_p->end()) {
    cache = &*q++;
  } else if (!use_cache) {
    cache = new FastMatCache;
  } else {
    branch_p->push_back(FastMatCache());
    cache = &branch_p->back();
  }
  return cache;
}
