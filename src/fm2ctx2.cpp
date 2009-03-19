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
::~CacheCtx2() { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx2
::Branch::Branch() : indx(-1) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2::jump(Branch &b) {
  if (b.indx==-1) {
    branchv.push_back(new clist_t());
    b.indx = branchv.size()-1;
  }
  branch_indx = b.indx;
  branch_p = branchv[branch_indx];
  q = branch_p->begin();
  was_cached = (q != branch_p->end());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::clear() { 
  for (unsigned int j=0; j<branchv.size(); j++) 
    delete branchv[j];
  branchv.clear(); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCache * 
FastMat2::CacheCtx2::step() {
  FastMatCache *cache=NULL;
  if (use_cache) {
    if (was_cached) {
      cache = &*q++;
      // was_cached = (q != branch_p->end());
    } else {
      branch_p->push_back(FastMatCache());
      cache = &branch_p->back();
      // was_cached = 0;
    }
//     printf("was_cached %d, cache %p, branch indx %d (%p)\n",
//            was_cached,cache,branch_indx,branch_p);
  } else {
    cache = new FastMatCache;
  }
  return cache;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2::print() {
  for (unsigned int j=0; j<branchv.size(); j++) {
    printf("branch indx %d, ptr %p: ",j,branchv[j]);
    clist_t &l = *branchv[j];
    clist_t::iterator q1 = l.begin();
    while (q1 != l.end()) printf(" %p",&*q1++);
    printf("\n");
  }
}
