///__INSERT_LICENSE__
//$Id mstorti-v6-branch-1.0.1-5-ge86f38c Wed Sep 19 13:06:15 2007 -0300$

#include <cmath>
#include <cstdio>

using namespace std;
#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/autostr.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
#if 0
FastMatCache *
FastMat2::CacheCtx::step(const char *label,
                         const FastMat2 *p1,
                         const FastMat2 *p2,
                         const FastMat2 *p3) { 
  return step(); 
}
#endif

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
void FastMat2::CacheCtx2
::check_clear() { as.clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::check(const char *label) { 
  as.sprintf("%s ",label);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::check(const FastMat2 *Ap) { 
  as.sprintf("%p ",Ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::check(const char *label,const FastMat2 *Ap,
        const FastMat2 *Bp,const FastMat2 *Cp) {
  check(label);
  if (Ap) check(Ap);
  if (Bp) check(Bp);
  if (Cp) check(Cp);
}

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
FastMatCache* FastMat2::CacheCtx2
::step() {
  FastMatCache *cache=NULL;
  if (use_cache) {
    if (was_cached) {
      cache = &*q++;
      if (do_check_labels)
        PETSCFEM_ASSERT(as.str()==cache->check_label,
                        "Failed FastMat2 cache check, "
                        "cached: \"%s\", wanted: \"%s\"",
                        cache->check_label.c_str(),as.str()); 
    } else {
      branch_p->push_back(FastMatCache());
      cache = &branch_p->back();
      if (do_check_labels)
        cache->check_label = as.str();
    }
    as.clear();
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
::check_labels(int do_check) {
  do_check_labels = do_check; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx2
::Branchv::init(int d0,int d1,int d2,int d3) {
  PETSCFEM_ASSERT0(shape.size()==0,"cant't resize a branch vector");  
  PETSCFEM_ASSERT0(rank==0,"cant't resize a branch vector");  
  if (d0>0) shape.push_back(d0);
  if (d1>0) shape.push_back(d1);
  if (d2>0) shape.push_back(d2);
  if (d3>0) shape.push_back(d3);
  rank = shape.size();
  unsigned int sz = 1;
  for (int j=0; j<rank; j++) sz *= shape[j];

  PETSCFEM_ASSERT0(sz>0,"total number of branchs in "
                   "branch vector must be positive");  
  
  bv.resize(sz);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx2
::Branchv::Branchv(int d0,int d1,int d2,int d3) 
  : rank(0) { if (d0>=0) init(d0,d1,d2,d3); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx2
::Branchv::Branchv(const vector<int> &shape_a) {
  shape = shape_a;
  rank = bv.size();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx2::Branch &
FastMat2::CacheCtx2::Branchv
::operator()(int j0,int j1,int j2,int j3 ) {
  int indx=j0;
  if (j1>=0) indx = indx*shape[1]+j1;
  if (j2>=0) indx = indx*shape[2]+j2;
  if (j3>=0) indx = indx*shape[3]+j3;
  PETSCFEM_ASSERT0(indx>=0,"indx");
  return bv[indx];
}
