//__INSERT_LICENSE__
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <unistd.h>

#include <map>
#include <algorithm>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fastlib2.h>
#include <src/fm2prod.h>
#include <src/vecmacros.h>

int FASTMAT2_USE_PROD2=1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::get_addresses(Indx permA,Indx Afdims,
                             vector<double *> &ap) const {
  int n = permA.size();
  if (ap.empty()) return;
  Indx iA(n,1), aindx(n,1),dA(n,-1);
  for (int j=0; j<n; j++) dA[j] = Afdims[permA[j]];
  int j=0;
  while (1) {
    for (int l=0; l<n; l++) {
      int k = permA[l];
      aindx[k] = iA[l];
    }
    ap[j++] = location(aindx);
    if (!inc(iA,dA)) break;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static void 
check_superlinear(vector<double *> &ap, int nrow,int ncol,
                  int &ok,int &lda,
                  CBLAS_TRANSPOSE &trans,int trans_ok) {
  int inccol = 1;
  int incrow = ncol;
  int sz = nrow*ncol;
  if (sz && ncol>1) inccol = int(ap[1]-ap[0]);
  if (sz && nrow>1) incrow = int(ap[ncol]-ap[0]);
  ok=1;
  for (int j=0; j<nrow; j++) {
    for (int k=0; k<ncol; k++) {
      if (int(ap[j*ncol+k]-ap[0]) != j*incrow+k*inccol) {
        ok=0;
        break;
      }
    }
  }
  lda = ncol;
  trans = CblasNoTrans;
  if (!ok) return;
  if (inccol==1) {
    trans = CblasNoTrans;
    lda = incrow;
    return;
  } 
  if (trans_ok && incrow==1) {
    trans = CblasTrans;
    lda = inccol;
    return;
  }
  ok = 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod2_subcache_t::init() {
  if (!table_was_initialized) {
    table_was_initialized = 1;
    load_funs();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int prod2_subcache_t
::table_indx(int n,int m,int p,int jat,int jbt) {
  return VEC_ADDR_5(n-1,m-1,NMAX,p-1,NMAX,jat,2,jbt,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
prod2_subcache_t::gemm_fun_t prod2_subcache_t
::get_fun(int n,int m,int p,int jat,int jbt) {
  int indx = table_indx(n,m,p,jat,jbt);
  return table[indx];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod2_subcache_t
::table_load(int n,int m,int p,int jat,int jbt,gemm_fun_t f) {
  table[table_indx(n,m,p,jat,jbt)] = f;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod2_subcache_t
::init(const FastMat2 &A,const FastMat2 &B,FastMat2 &C,
       vector<int> &ixa,vector<int> &ixb) {
  init();
  // get the free indices of A and B
  Indx ii,Afdims,Bfdims,Cfdims;
  // Afdims, Bfdims:= are the dimensions of the `free'
  // indices of A,B (not free indices of the contraction,
  // but the indices of A,B that are not restricted)
  A.get_dims(Afdims);
  B.get_dims(Bfdims);
  
  // Nbr of indices in A and B
  int niA = Afdims.size();
  int niB = Bfdims.size();
  // Total number of input indices
  int ndims = niA+niB;

  // Check the indices that were passed for
  // the contraction are OK in length
  int ixas = int(ixa.size());
  int ixbs = int(ixb.size());
  PETSCFEM_ASSERT(ixas == niA,
                  "Nbr of indices passed in the contraction for A"
                  " must be equal to the actual indices of A. "
                  "nbr of A indices %d. Passed %d",ixas,niA);  
  PETSCFEM_ASSERT(ixbs == niB,
                  "Nbr of indices passed in the contraction for B"
                  " must be equal to the actual indices of B. "
                  "nbr of B indices %d. Passed %d",ixbs,niB);  
    
  // This implementation is based on an older one
  // that used the vector of indices for contraction.
  // Then we simply put the `ixa' and `ixb' indices in `ii'
  // and continue. 
  for (unsigned int j=0; j<ixa.size(); j++) 
    ii.push_back(ixa[j]);
  for (unsigned int j=0; j<ixb.size(); j++) 
    ii.push_back(ixb[j]);

  // If `C' is not defined we have to initialize it
  if (!C.defined) {
    // Count total number of free (not contracted) indices
    int nf=0;
    for (int j=0; j<ndims; j++) nf += ii[j]>0;
    // Dimension index vector
    Indx ndimsf(nf,0);
    // Get dimensions
    for (int j=0; j<ndims; j++) {
      int k = ii[j];
      if (k>0) {
        PETSCFEM_ASSERT(k<=nf,
                        "Index at position %d is %d, but "
                        "nbr of positive (free) indices is smaller: %d\n",
                        j,k,nf);  
        ndimsf[k-1] = (j<niA? Afdims[j] : Bfdims[j-niA]);
      }
    }
    C.create_from_indx(ndimsf);
  }

  // Get free dims of output matrix
  C.get_dims(Cfdims);
  int niC = Cfdims.size();

  // Nbr of free and contracted indices for A and B
  int nfA=0,nfB=0,nfree=0,nc=0;
  for (int j=0; j<niA; j++) {
    if (ixa[j]>0) nfA++;
    else nc++;
  }
  for (int j=0; j<niB; j++) {
    if (ixb[j]>0) nfB++;
    else nc++;
  }
  PETSCFEM_ASSERT(nc%2==0,
                  "Nbr of contracted indices must be even."
                  "nc  %d",nc);  

  nc /= 2;
  nfree = nfA+nfB;
  int ntot = niA+niB+niC;
  PETSCFEM_ASSERT(ntot==2*(nfree+nc),
                  "Bad balance of indices. "
                  "niA %d, niB %d, niC %d, nfree %d, nctr %d",
                  niA,niB,niC,nfree,nc);  
  

  // We reorder the indices of A in (fA,cA), B in (cB,fB). 
  // `f' denotes free and `c' contracted.
  // Indices in fA and fB keep the order as was given.
  // cA and cB are reordered according to the negative
  // index that was given. The indices in C are also
  // reordered so as to keep the same order as appear
  // in fA,fB concatenated.
  // Example: if we make C(lim)=A(kij)B(lkmj) 
  // then: we reorder this as
  // C(ilm)=A(ijk)B(kjlm) 
  // So that we have the following permutations
  // mapC=[1,0,2] mapA=[2,0,1]
  // At the same time we compute the
  // row dimension of A = product of free indices,
  // and column dimension = product of contracted indices. 
  // This is as if we convert A in a rectangular matrix
  // where all the free indices are colapsed in a row index
  // and the contracted indices in a single column index.
  // The same is done for B, but in this case they are inverted
  // contracted -> row, free -> column. 
  nrowa=1; ncola=1; nrowb=1; ncolb=1;
  Indx mapA(niA,-1),mapB(niB,-1),mapC(niC,-1);
  int kf=0;
  for (int j=0; j<niA; j++) {
    int k=ixa[j];
    if (k>0) {
      mapA[kf] = j;
      mapC[kf] = k-1;
      kf++;
      nrowa *= Afdims[j];
    } else {
      mapA[nfA-1-k] = j;
      ncola *= Afdims[j];
    }
  }
  nrowc = nrowa;
    
  kf=0;
  for (int j=0; j<niB; j++) {
    int k=ixb[j];
    if (k>0) {
      mapB[nc+kf] = j;
      mapC[nfA+kf] = k-1;
      kf++;
      ncolb *= Bfdims[j];
    } else {
      mapB[-1-k] = j;
      nrowb *= Bfdims[j];
    }
  }
  ncolc = ncolb;
  assert(ncola==nrowb);

  // Get dimensions of matrices
  nA = A.size();
  nB = B.size();
  nC = C.size();

  // Resize array of pointers.
  // ap is an array with the same size as A. 
  // Its elements are pointers to the elements
  // in A reordered as in the hypothetical
  // rectangular matrix. Same for B and C
  ap.resize(nA);
  bp.resize(nB);
  cp.resize(nC);

  // Get the addresses and put them in ap,bp,cp
  A.get_addresses(mapA,Afdims,ap);
  B.get_addresses(mapB,Bfdims,bp);
  C.get_addresses(mapC,Cfdims,cp);

  // Check if any of them are `superlinear'. This means
  // that elements in each row are equispaced by `inccol'
  // and the start of rows are equispaced by `incrow'
  check_superlinear(ap,nrowa,ncola,asl_ok,lda,transa,1);
  check_superlinear(bp,nrowb,ncolb,bsl_ok,ldb,transb,1);
  check_superlinear(cp,nrowa,ncolb,csl_ok,ldc,transc,1);

  // If `A' is not superlinear, then
  // its values are copied to the auxiliary matrices `a'
  // and we call dgemm on them. Same for the others.
  // So a,b,c must be allocated only if they are
  // not superlinear. 
  if (!asl_ok) a.resize(nA,0.0);
  if (!bsl_ok) b.resize(nB,0.0);
  if (!csl_ok) c.resize(nC,0.0);

  // We store in the cache a pointer `Ap' to either
  // the internal buffer (in fact to the first
  // element in the submatrix) OR to the
  // auxiliary buffer. 
  Ap=NULL; Bp=NULL; Cp=NULL; 
  if (nC>0) {
    Ap = (asl_ok? ap[0] : &a[0]);
    Bp = (bsl_ok? bp[0] : &b[0]);
    Cp = (csl_ok? cp[0] : &c[0]);
  }

  // Check for special fast functions
  use_fmgemm = 0;
  // nmax = 0; // Short-circuited to false!!
  int nmax1 = (NMAX<nmax ? NMAX : nmax);
  // if (!(nrowa>nmax1 || ncola>nmax1 || ncolb>nmax1 ||
  //       nrowa<1 || ncola<1 || ncolb<1 ||
  //       lda!=ncola || ldb!=ncolb || ldc!=ncolc || transc != CblasNoTrans)) {

  // Check if matrix dimensions are in range supported by FMGEMM
  int ok = nrowa<=nmax1 && ncola<=nmax1 && ncolb<=nmax1 &&
    nrowa>=1 && ncola>=1 && ncolb>=1;
  int jat = (transa==CblasNoTrans ? 0 : 1);
  int jbt = (transb==CblasNoTrans ? 0 : 1);
  int jct = (transc==CblasNoTrans ? 0 : 1);

  // Check appropriate leading dimensions in A and B
  ok = ok && (jat? lda==nrowa : lda==ncola);
  ok = ok && (jbt? ldb==nrowb : ldb==ncolb);
  ok = ok && (jct? ldc==nrowc : ldc==ncolc);
  if (ok) {
    use_fmgemm = 1;
    if (!jct) 
      gfun = table[table_indx(nrowa,ncola,ncolb,jat,jbt)];
    else
      gfun = table[table_indx(ncolb,ncola,nrowa,!jbt,!jat)];
  }
#ifdef DO_SIZE_STATS
  last_call_used_fmgemm = use_fmgemm;
  if (do_size_stats) {
    total_calls++;
    fmgemm_calls += use_fmgemm;
  }
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod2(const FastMat2 &A,const FastMat2 &B,
                vector<int> &ixa, 
                vector<int> &ixb) {
  return prod(A,B,ixa,ixb);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This is the function that makes the product of two matrices.
// The others for 3,4,etc... are wrappers to this one. 
// Computes C = A*B, where C is *this
FastMat2 & 
FastMat2::prod(const FastMat2 &A,const FastMat2 &B,
               vector<int> &ixa, 
               vector<int> &ixb) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod2",this,&A,&B);
    Indx Ixa,Ixb;
    for (unsigned int j=0; j<ixa.size(); j++)
      Ixa.push_back(ixa[j]);
    for (unsigned int j=0; j<ixb.size(); j++)
      Ixb.push_back(ixb[j]);
    ctx->check(Ixa);
    ctx->check(Ixb);
  }
#endif
  FastMatCache *cache = ctx->step();

  prod2_subcache_t *psc = NULL;
  if (!ctx->was_cached  ) {
    psc = new prod2_subcache_t;
    assert(psc);
    assert(!cache->sc);
    cache->sc = psc;
    // Make the initialization of this
    // product call
    psc->init(A,B,*this,ixa,ixb);
  }

  psc = dynamic_cast<prod2_subcache_t *>(cache->sc);
  assert(psc);
  psc->make_prod();
  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static CBLAS_TRANSPOSE flip(CBLAS_TRANSPOSE flag) {
  return (flag==CblasNoTrans? CblasTrans : CblasNoTrans);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
const char *cbt2str(CBLAS_TRANSPOSE flag) {
  return (flag==CblasTrans? "T" : "N");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod2_subcache_t::make_prod() {
  // Makes the product by calling to BLAS DGEMM
  if (nC==0) return;

#ifdef DO_SIZE_STATS
  double start=NAN;
  if (do_size_stats) start = MPI_Wtime();
#endif

  // Copy from A,B internal buffers to
  // auxiliary buffers if needed
  if (!asl_ok) for (int j=0; j<nA; j++) a[j] = *ap[j];
  if (!bsl_ok) for (int j=0; j<nB; j++) b[j] = *bp[j];

  // printf("use FMGEMM: %d\n",FASTMAT2_USE_FMGEMM && use_fmgemm);
  if (FASTMAT2_USE_FMGEMM && use_fmgemm) {
    if (transc==CblasNoTrans) gfun(Ap,Bp,Cp);
    else gfun(Bp,Ap,Cp);
#if 0
    // Perform a check
    vector<double> cv(nrowc*ncolc);
    if (transc==CblasNoTrans) 
      cblas_dgemm(CblasRowMajor,transa,transb,nrowa,ncolb,
                  ncola,1.0,Ap,lda,Bp,ldb,0.0,&cv[0],ldc);
    else cblas_dgemm(CblasRowMajor,flip(transb),flip(transa),
                     ncolb,nrowa,ncola,1.0,Bp,ldb,Ap,lda,0.0,
                     &cv[0],ldc);
    double erro=0.0,tol=1e-10,cnorm=0.0;
    for (unsigned int j=0; j<cv.size(); j++) {
      erro += fabs(cv[j]-Cp[j]);
      cnorm += Cp[j];
    }
    printf("FMGEMM error %g, tol %g, ||C|| %g\n",erro,tol,cnorm);
    // PETSCFEM_ASSERT(erro<tol,"FMGEMM error %g, tol %g, ||C|| %g\n",erro,tol,cnorm);
#endif
  } else {
    // Call DGEMM
#ifdef CALL_BLAS_LAPACK_FROM_PETSC
    PetscScalar _DOne=1.0,_DZero=0.0;

    // This way it works, but makes C'=A'*B'
    if (transc==CblasNoTrans) 
      BLASgemm_(cbt2str(transa),cbt2str(transb),
                &nrowa,&ncolb,&ncola,
                &_DOne,Ap,&lda,
                Bp,&ldb,
                &_DZero,Cp,&ldc);
    else 
      BLASgemm_(cbt2str(flip(transb)),cbt2str(flip(transa)),
                &ncolb,&nrowa,&ncola,
                &_DOne,Bp,&ldb,
                Ap,&lda,
                &_DZero,Cp,&ldc);
#else
    if (transc==CblasNoTrans) 
      cblas_dgemm(CblasRowMajor,transa,transb,nrowa,ncolb,
                  ncola,1.0,Ap,lda,Bp,ldb,0.0,Cp,ldc);
    else cblas_dgemm(CblasRowMajor,flip(transb),flip(transa),
                     ncolb,nrowa,ncola,1.0,Bp,ldb,Ap,lda,0.0,
                     Cp,ldc);
#endif
  }

  // Copy from c auxiliary buffer to C
  // internal buffer if needed
  if (!csl_ok) for (int j=0; j<nC; j++) *cp[j] = c[j];

#ifdef DO_SIZE_STATS
  if (do_size_stats) {
    mat_sz_t m;
    m.m = nrowa; m.n=ncola; m.p=ncolb;
    stats_t &s = stat_table[m];
    s.count++;
    s.time += MPI_Wtime()-start;
  }
#endif

}

#include "./fmgemmdefs.h"  
int prod2_subcache_t::NMAX=PF_FMGEMM_NMAX;
int prod2_subcache_t::nmax=PF_FMGEMM_NMAX;
int prod2_subcache_t::table_was_initialized=0;
int prod2_subcache_t::FASTMAT2_USE_FMGEMM=1;
vector<prod2_subcache_t::gemm_fun_t> prod2_subcache_t::table;

typedef void (*gemm_fun_t)(double *a,double *b,double *c);
extern void prod2_subcache_load_funs(vector<gemm_fun_t> &table);

void prod2_subcache_t::load_funs() {
  table.resize(4*NMAX*NMAX*NMAX);
  prod2_subcache_load_funs(table);
}

#ifdef DO_SIZE_STATS

bool prod2_subcache_t
::comp(const prod2_subcache_t::spair_t &a,
       const prod2_subcache_t::spair_t &b) {
  return a.second.time>b.second.time;
}

int prod2_subcache_t::total_calls=0;
int prod2_subcache_t::fmgemm_calls=0;
int prod2_subcache_t::do_size_stats=0;
int prod2_subcache_t::last_call_used_fmgemm=-1;

map<prod2_subcache_t::mat_sz_t,prod2_subcache_t::stats_t> prod2_subcache_t::stat_table;
void prod2_subcache_t::report_stats() {
  map<mat_sz_t,stats_t>::iterator q = stat_table.begin();
  double ctime=0.0;
  vector<spair_t> v;
  while (q!=stat_table.end()) {
    v.push_back(spair_t(q->first,q->second));
    q++;
  }
  sort(v.begin(),v.end(),comp);
  double total = 0.0;
  for (unsigned int j=0; j<v.size(); j++) {
    const stats_t &s = v[j].second;
    total += s.time;
  }
  for (unsigned int j=0; j<v.size(); j++) {
    const mat_sz_t &z = v[j].first;
    const stats_t &s = v[j].second;
    ctime += s.time;
    printf("%d %d %d %d %g [%.2f%%] (cum %g [%f%%])\n",
           z.m,z.n,z.p,s.count,
           s.time,100.0*s.time/total,
           ctime,ctime/total*100.0);
  }
  printf("Total prods %d, use FMGEMM %d\n",total_calls,fmgemm_calls);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
bool prod2_subcache_t::mat_sz_t
::operator<(const mat_sz_t &rhs) const {
  if (m!=rhs.m) return m<rhs.m;
  if (n!=rhs.n) return n<rhs.n;
  return p<rhs.p;
}
#endif
