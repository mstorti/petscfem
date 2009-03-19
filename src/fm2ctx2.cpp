///__INSERT_LICENSE__
//$Id mstorti-v6-branch-1.0.1-5-ge86f38c Wed Sep 19 13:06:15 2007 -0300$

#include <cmath>
#include <cstdio>

using namespace std;
#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/autostr.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCache *
FastMat2::CacheCtx::step(const char *label,
                         const FastMat2 *p1,
                         const FastMat2 *p2,
                         const FastMat2 *p3) { 
  return step(); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2
::CacheCtx2::CacheCtx2() :
  branch_indx(-1), 
  branch_p(NULL),
  do_check_labels(0) { }

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
  return step(NULL);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCache* FastMat2::CacheCtx2
::step(const char *label, const FastMat2 *p1,
       const FastMat2 *p2, const FastMat2 *p3) {
  FastMatCache *cache=NULL;
  if (use_cache) {
    AutoString as;
    if (do_check_labels) {
      as.sprintf("%s %p %p %p",label,p1,p2,p3);
      printf("check_label: %s\n",as.str());
    }
    
    if (was_cached) {
      cache = &*q++;
      if (do_check_labels)
        PETSCFEM_ASSERT(as.str()==cache->check_label,
                        "Failed FastMat2 cache check"
                        "Cached \"%s\", wanted \"%s\"",
                        cache->check_label.c_str(),as.str()); 
    } else {
      branch_p->push_back(FastMatCache());
      cache = &branch_p->back();
      cache->check_label = as.str();
    }
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::check_labels() { do_check_labels = 1; }
