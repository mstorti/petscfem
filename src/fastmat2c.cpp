//__INSERT_LICENSE__
#include <cmath>
#include <ctime>
#include <cstdio>
#include <unistd.h>

#include <algorithm>
#ifdef USE_MKL
#include <mkl_cblas.h>
#endif

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fastlib2.h>
#include <src/fm2prod.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
mat_info::mat_info() 
  : Ap(NULL), 
    type(multiprod_subcache_t::UNKNOWN),
    is_active(multiprod_subcache_t::UNDEF) {}

// The definitions are generated with a perl
// script `mprod.pl' in a header and included here. 
#include "./mproddef.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Searches for `x' in vector `v'
static int vfind(const vector<int> &v,int x) {
  int n=v.size();
  for (int j=0; j<n; j++)
    if (v[j]==x) return 1;
  return 0;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Auxiliary matrix to print th information for the product
static void 
print_mat_info(mat_info_cont_t::iterator q,
               const char *s=NULL) {
  if (s) printf("%s: ",s);
  mat_info &mi = q->second;
  vector<int> 
    &qc = mi.contract,
    &qd = mi.dims;
  printf("key: %d, ptr %p, type %d, dims: ",
         q->first,mi.Ap,mi.type);
  int rank = qd.size();
  for (int j=0; j<rank; j++) printf("%d ",qd[j]);
  printf("ctr-pos: ");
  for (int j=0; j<rank; j++) printf("%d ",qc[j]);
  printf("\n");
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static string print_label(const string &label) {
  int n= label.size();
  assert(n>=2);
  return label.substr(1,n-2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Computes and returns the operation count for a
// possible product. Also returns the product of the free
// indices in q and r, and of the contracted indices.
// (BTW the opcount is qfree*rfree*
int compute_opcount(const mat_info &qmi,
                    const mat_info &rmi,
                    mat_info &smi,
                    int &qfree,int &rfree,int &qctr) {
  const vector<int> 
    &qc = qmi.contract,
    &qd = qmi.dims,
    &rc = rmi.contract,
    &rd = rmi.dims;
  vector<int> 
    &sc = smi.contract,
    &sd = smi.dims;
  int 
    rctr=1,
    qrank = qd.size(),
    rrank = rd.size(),
    k,dim;
  rfree=1;
  qfree=1;
  qctr=1;
  for (int j=0; j<qrank; j++) {
    k = qc[j];
    dim = qd[j];
    // The free indices are those that
    // are positive or negative but they are
    // not in the other matrix. There may be negative
    // indices because they are contracted in other
    // product.
    if (k>0 || !vfind(rc,k)) { 
      qfree *= dim;
      sc.push_back(k);
      sd.push_back(dim);
    } else qctr *= dim;
  }
  for (int j=0; j<rrank; j++) {
    k = rc[j];
    dim = rd[j];
    // Same as before
    if (k>0 || !vfind(qc,k)) {
      rfree *= dim;
      sc.push_back(k);
      sd.push_back(dim);
    }
    else rctr *= dim;
  }
  assert(qctr==rctr);
  return qfree*rfree*qctr;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
multiprod_subcache_t::~multiprod_subcache_t() {
  // Deallocates the temporary matrices
  int n = mat_info_cont.size();
  for (int j=0; j<n; j++) {
    mat_info &mi = mat_info_cont[j];
    if (mi.type == TMP && mi.Ap) delete mi.Ap;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void multiprod_subcache_t::make_prod() {
  // The `order' table contained for each product
  // 3 integers: the destination matrix, and the 2
  // input matrix
  int nprod = nmat-1,
    qkey, rkey, skey;
    //#define VERBOSE
#ifdef VERBOSE
  char line[100];
#endif
  for (int j=0; j<nprod; j++) {
    mat_info 
      // Destination matrices
      &smi = mat_info_cont[order[3*j]],
      // Source matrices
      &qmi = mat_info_cont[order[3*j+1]],
      &rmi = mat_info_cont[order[3*j+2]];

    // Makes the product, calling prod(a,b,indxa,indxb)
    smi.Ap->prod(*qmi.Ap,*rmi.Ap,qmi.contract,rmi.contract);

#ifdef VERBOSE
    sprintf(line,"q[%d]: ",order[3*j+1]);
    qmi.Ap->print(line);

    sprintf(line,"r[%d]: ",order[3*j+2]);
    rmi.Ap->print(line);

    sprintf(line,"s[%d]: ",order[3*j]);
    smi.Ap->print(line);
#endif
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// General case. Receives a list of pointers to matrices
// and a vector of indices. The instances for 2,3,4,... matrices
// are wrappers to this. In turn this routine only computes the
// best order to do the products and calls the version
// for two matrices
FastMat2 & 
FastMat2::prod(vector<const FastMat2 *> &mat_list,
               Indx &indx) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_indx",this);
    int nmat = mat_list.size();
    for (int j=0; j<nmat; j++)
      ctx->check(mat_list[j]);
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  multiprod_subcache_t *mpsc=NULL;
  if (!ctx->was_cached  ) {
    mpsc = new multiprod_subcache_t(cache,ctx,*this,mat_list,indx);
    assert(!cache->sc);
    cache->sc = mpsc;
  }

  mpsc = dynamic_cast<multiprod_subcache_t *> (cache->sc);
  assert(mpsc);
  mpsc->make_prod();

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
multiprod_subcache_t
::multiprod_subcache_t(FastMatCache *cache_a,FastMat2::CacheCtx *ctx,
                       FastMat2 &result,
                       vector<const FastMat2 *> &mat_list,
                       Indx &indx) {
  // Check sum of ranks of matrices equals
  // size of index vector passed
  nmat = mat_list.size();
  PETSCFEM_ASSERT0(nmat>=2,"prod needs at least 2 matrices");  

  int nindx = 0;
  for (int j=0; j<nmat; j++)
    nindx += mat_list[j]->n();
  PETSCFEM_ASSERT(nindx==indx.size(),"Not correct nbr of indices. "
                  "Sum of ranks %d. nbr of indices %d",nindx,indx.size());  

  // For contracted indices (negative) count how many times
  // they appear (should be 2).
  // For the free indices they should appear once and
  // in the range [1,rank] where rank is the
  // rank of the output matrix.
  map<int,int> indices; 
  for (int j=0; j<nindx; j++) {
    indices[indx[j]]++;
  }

  map<int,int>::iterator s = indices.begin();
  int erro=0, outrank=-1;
  while (s!=indices.end()) {
    int k=s->first, n=s->second;
    if (k<0 && n!=2) {
      erro=1;
      printf("contracted index %d repeated %d times (should be 2)\n",k,n);
    } 
    if (k>0 && k > outrank) outrank = k;
    if (k>0 && n!=1) {
      printf("free index %d repeated %d times (should be 1)\n",k,n);
      erro = 1;
    }
    if (k==0) printf("index null is invalid %d\n",k);
    s++;
  }
  
  // Check all free indices are in range [1,rank]
  int erro1=0;
  for (int j=1; j<=outrank; j++) {
    if (indices.find(j)==indices.end()) erro1=1;
  }
  if (erro1) {
    erro = 1;
    printf("outrank %d and following indices are missing: ",outrank);
    for (int j=1; j<=outrank; j++) 
      if (indices.find(j)==indices.end()) printf("%d ",j);
    printf("\n");
  }
  if (erro) PETSCFEM_ERROR0("detected errors, aborting\n");  

  // Stores the info of matrices (contracted indices and
  // pointer to matrices) in a list of structures mat_info
  // vector<mat_info> &mat_info_cont = mpsc->mat_info_cont;
  mat_info_cont.resize(nmat);
  // vector<int> &order = mpsc->order;
  PETSCFEM_ASSERT0(order.empty(),"Internal error");  
     
  // Load the input matrices in the container
  int j=0;
  for (int k=0; k<nmat; k++) {

    // Pointer to the matrix
    FastMat2 &A = *(FastMat2 *)mat_list[k];
    mat_info &m = mat_info_cont[k];
    m.Ap = &A;
    m.type = multiprod_subcache_t::OLD;
    m.position = k;
    int rank = A.n();
    // This stores the indices to be contracted
    m.contract.resize(rank);
    // This stores the dims
    m.dims.resize(rank);
    for (int l=0; l<rank; l++) {
      m.contract[l] = indx[j++];
      m.dims[l] = A.dim(l+1);
    }
  }

  // In `mat_info_cont_plain' the mat_infos are stored
  // in compact form, i.e. when matrix in positions `q'
  // and `r' are contracted `q' is replaced with the product
  // and `r' is removed
  mat_info_cont_t 
    mat_info_cont_plain = mat_info_cont;
  vector<string> labels(nmat);
  if (fastmat_stats.labels) {
    for (int j=0; j<nmat; j++) {
      assert(fastmat_stats.labels[j]);
      labels[j] = fastmat_stats.labels[j];
    }
    fastmat_stats.nmat = nmat;
    assert(!fastmat_stats.labels[nmat]);
  }
  vector<int> opt_order;
  intmax_t nopso=-1;
  // double start=MPI_Wtime();
  if (ctx->mprod_order==FastMat2::CacheCtx::optimal
      || (ctx->mprod_order==FastMat2::CacheCtx::mixed 
          && nmat <= ctx->optimal_mprod_order_max)) {
    nopso = compute_optimal_order(mat_info_cont_plain,opt_order);
  } else if (ctx->mprod_order==FastMat2::CacheCtx::heuristic
             || (ctx->mprod_order==FastMat2::CacheCtx::mixed 
                 && nmat > ctx->optimal_mprod_order_max)) {
    nopso = compute_heuristic_order(mat_info_cont_plain,opt_order);
  } else if (ctx->mprod_order==FastMat2::CacheCtx::natural 
             || ctx->mprod_order==FastMat2::CacheCtx::reverse) {
    int reverse = ctx->mprod_order==FastMat2::CacheCtx::reverse;
    nopso = compute_natural_order(mat_info_cont_plain,
                                  opt_order,reverse);
  } else {
    PETSCFEM_ERROR("unknown mprod_order %d",ctx->mprod_order);  
  }
  // printf("nmat %d, elapsed %g secs\n",nmat,MPI_Wtime()-start);
  assert(int(opt_order.size())==2*(nmat-1));
  ctx->mprod_nopscount = nopso;

#if 0
  for (int j=0; j<nmat-1; j++) {
    char label[1000];
    int q=opt_order[2*j], r=opt_order[2*j+1];
    sprintf(label,"(%s*%s)",labels[q].c_str(),labels[r].c_str());
    labels[q] = label;
    labels.erase(labels.begin()+r);
  }
  printf("OPTIMAL product order: %s\n",print_label(labels[0].c_str()));
#endif

#if 0
  for (int j=0; j<nmat-1; j++) {
    printf("prod %d: a%d = a%d * a%d\n",
           j,opt_order[2*j],opt_order[2*j],opt_order[2*j+1]);
  }
  exit(0);
#endif

  // The created temporaries are assigned
  // this key (position in `mat_info_cont')
  int new_mat_indx = nmat;
  // These are the keys stored in the compact version
  vector<int> keys(nmat);
  for (int j=0; j<nmat; j++) keys[j] = j;

  int q,r,qkey,rkey,qrank, rrank, k, dim;
  // Start simulating the products and storing
  // the info in `order'
  for (int j=0; j<nmat-1; j++) {
    int nact = nmat-j;

    // Insert a new `mat_info' corresponding to
    // the temporary to be created
    mat_info_cont.push_back(mat_info());
    int skey = new_mat_indx++;
    mat_info &smi = mat_info_cont[skey];
    // For the firs products except the last one
    // we reorder the free indices
    // so that they are in range [1,n].
    // For the last one this is not needed because
    // it should end in this way automatically.
    // And the order should be preserved becuase
    // otherwise the elements in the resulting
    // matrix would be in an order not desired. 
    int reo=1;
    if (nact>2) {
      // This is for all the products except the last one
      smi.Ap = new FastMat2(ctx);
      smi.type = multiprod_subcache_t::TMP;
    } else {
      // In this case the output matrix is *this
      smi.Ap = &result;
      // It is NOT a temporary
      smi.type = multiprod_subcache_t::OLD;
      // Do not reorder the indices
      reo = 0;
    }
    // This new matrix is active
    smi.is_active = multiprod_subcache_t::ACTIVE;
    vector<int> 
      &sc = smi.contract,
      &sd = smi.dims;

    q = opt_order[2*j];
    r = opt_order[2*j+1];
    assert(q<r);
    qkey = keys[q];
    rkey = keys[r];

    // Product to be done is S = Q * R
    // get info for Q matrix
    mat_info &qmi = mat_info_cont[qkey];
    qmi.is_active = multiprod_subcache_t::INACTIVE;
    vector<int> 
      &qc = qmi.contract,
      &qd = qmi.dims;

    // get info for R matrix
    mat_info &rmi = mat_info_cont[rkey];
    rmi.is_active = multiprod_subcache_t::INACTIVE;
    vector<int> 
      &rc = rmi.contract,
      &rd = rmi.dims;

    qrank = qd.size();
    rrank = rd.size();
    int sindx = 1;
    // Compute the new contraction indices for Q
    for (int j=0; j<qrank; j++) {
      k = qc[j];
      dim = qd[j];
      if (k>0 || !vfind(rc,k)) {
        if (reo) qc[j] = sindx++;
        sc.push_back(k);
        sd.push_back(dim);
      }
    }
    // Compute the new contraction indices for R
    for (int j=0; j<rrank; j++) {
      k = rc[j];
      dim = rd[j];
      if (k>0 || !vfind(qc,k)) {
        if (reo) rc[j] = sindx++;
        sc.push_back(k);
        sd.push_back(dim);
      }
    }
    // Insert the `keys' for this product in `order'
    order.push_back(skey);
    smi.position = qmi.position;
    order.push_back(qkey);
    order.push_back(rkey);
    // Update `keys' vector
    keys[q] = skey;
    keys.erase(keys.begin()+r);
  }

#if 1
  if (fastmat_stats.labels 
      &&fastmat_stats.print_prod_order) {
    int sindx = nmat;
    vector<string> labels(nmat);
    for (int j=0; j<nmat; j++) {
      assert(fastmat_stats.labels[j]);
      labels[j] = fastmat_stats.labels[j];
    }
    for (int j=0; j<nmat-1; j++) {
      char label[1000];
      int q=order[3*j+1], r=order[3*j+2];
      sprintf(label,"(%s*%s)",labels[q].c_str(),labels[r].c_str());
      labels.push_back(label);
    }
    assert(int(labels.size())==2*nmat-1);
    printf("product order: %s\n",
           print_label(labels[2*nmat-2]).c_str());
  }
#endif

  assert(order.size() == (unsigned int)(3*(nmat-1)));

  // Now for each product to be done checks that
  // all negative indices are `packed', i.e.
  // in range [1,n] with all indices used.
  for (int j=0; j<nmat-1; j++) {
    mat_info 
      &qmi = mat_info_cont[order[3*j+1]],
      &rmi = mat_info_cont[order[3*j+2]];
    vector<int> 
      &qc = qmi.contract,
      &rc = rmi.contract;

    // new_neg maps the `old' contraction indices
    // to the `new' indices
    map<int,int> new_neg;
    // this is the next new negative indices to
    // be assigned
    int new_neg_free=-1;
    int k,
      qrank = qc.size(),
      rrank = rc.size();
    for (int l=0; l<qrank; l++) {
      k = qc[l];
      if (k<0) {
        // If another negative index appears
        // we assign a new index
        if (new_neg.find(k)==new_neg.end()) 
          new_neg[k] = new_neg_free--;
        // replace the old indices with the
        // new indices
        qc[l] = new_neg[k];
      }
    }
    // Same for the q indices
    for (int l=0; l<rrank; l++) {
      k = rc[l];
      if (k<0) {
        if (new_neg.find(k)==new_neg.end()) 
          new_neg[k] = new_neg_free--;
        rc[l] = new_neg[k];
      }
    }
#ifdef VERBOSE
    printf("prod #%d [qc:",j);
    for (int l=0; l<qrank; l++) 
      printf("%d,",qc[l]);
    printf("][rc:");
    for (int l=0; l<rrank; l++) 
      printf("%d,",rc[l]);
    printf("]\n");
#endif
  }
}

fastmat_stats_t fastmat_stats;

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This is the function that makes the product of two matrices.
// The others for 3,4,etc... are wrappers to this one. 
FastMat2 & 
FastMat2::prod(const FastMat2 &A,const FastMat2 &B,
               vector<int> &ixa, 
               vector<int> &ixb) {

#if 0
  static int reported=0;
  if (!reported) {
    printf("FASTMAT2_USE_PROD2=%d\n",FASTMAT2_USE_PROD2);
    const char *var = getenv("OMP_NUM_THREADS");
    printf("OMP_NUM_THREADS %s\n",(var? var: "unknown"));
    reported=1;
  }
#endif

  if (FASTMAT2_USE_PROD2) {
    prod2(A,B,ixa,ixb);
    return *this;
  }

  double start=MPI_Wtime();
  fastmat_stats.ncall++;

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

  prod_subcache_t *psc=NULL;
  if (!ctx->was_cached  ) {
    fastmat_stats.ncall_not_cached++;
    Indx ia,ib,ii;

    // get the free indices of A and B
    Indx Afdims,Bfdims,fdims;
    A.get_dims(Afdims);
    B.get_dims(Bfdims);

    int niA = Afdims.size();
    int niB = Bfdims.size();
    int ndims = niA+niB;
    
    // This implementation is based on an older one
    // that used the vector of indices for contraction.
    // Then we simply put the `ixa' and `ixb' indices in `ii'
    // and continue. 
    for (unsigned int j=0; j<ixa.size(); j++) 
      ii.push_back(ixa[j]);
    for (unsigned int j=0; j<ixb.size(); j++) 
      ii.push_back(ixb[j]);

    // Count how many free and contracted indices are
    int nfree=0,nc=0;
    for (int j=0; j<ii.size(); j++) {
      int k = ii[j];
      if (k>0 && k>nfree) nfree = k;
      if (k<0 && -k>nc) nc = -k;
    }

    // Get the free and contracted indices
    Indx ifree(nfree,0),icontr(2*nc,0);
    for (int j=0; j<ii.size(); j++) {
      int k = ii[j];
      if (k>0) {
	ifree[k-1] = j+1;
      } else {
	k = -k;
	if (icontr[2*(k-1)]==0) {
	  icontr[2*(k-1)]=j+1;
	} else {
	  icontr[2*k-1]=j+1;
	}
      }
    }
    // ifree.print("ifree: ");
    // icontr.print("icontr: ");
    
    // ndimsf:= dimensions of the free indices
    // ndimsc:= dimensions of the contracted indices
    Indx ndimsf(nfree,0),ndimsc(nc,0);
    for (int j=0; j<nfree; j++) {
      int k = ifree[j];
      if (k<=niA) {
	ndimsf[j] = Afdims[k-1];
      } else {
	ndimsf[j] = Bfdims[k-niA-1];
      }
    }

    // ndimsf.print("dimensions of the free part: ");

    // Dimension C if necessary
    if (!defined) create_from_indx(ndimsf);
    // If the resulting matrix is 2-size then do nothing
    if (comp_storage_size(ndimsf)==0) {
      cache->nelems=0;
      return *this;
    }

    // Get free dims of output matrix
    get_dims(fdims);
    if (ndimsf != fdims) {
      Afdims.print("A free dims: ");
      Bfdims.print("B free dims: ");
      ndimsf.print("Combined free dims: ");
      fdims.print("Free dims on result matrix: ");
      PETSCFEM_ERROR0("Combined free dims doesn't match free"
	     " dims of  result.\n");
    }

    // Check that the paired contracted indices are
    // equal.
    for (int j=0; j<nc; j++) {
      int k1 = icontr[2*j];
      int k2 = icontr[2*j+1];
      int nd1,nd2;
      if (k1<=niA) {
	nd1 = Afdims[k1-1];
      } else {
	nd1 = Bfdims[k1-niA];
      }
      if (k2<=niA) {
	nd2 = Afdims[k2-1];
      } else {
	nd2 = Bfdims[k2-niA-1];
      }
      PETSCFEM_ASSERT0(nd1==nd2,"paired contracted indices don't match");  
      ndimsc[j]=nd1;
    }

    // ndimsc.print("dimensions of the contracted part: ");

    Indx findx(nfree,1),cindx(nc,1),tot_indx(ndims,0),
      aindx(niA,0),bindx(niB,0);

    // Loading addresses in cache. For each element in the
    // distination target, we store the complete list of
    // addresses of the lines of elements that contribute to
    // it.
    cache->nlines = mem_size(ndimsf);
    cache->prod_cache.resize(cache->nlines);
    cache->line_cache_start = &*cache->prod_cache.begin();
    cache->line_size = mem_size(ndimsc);
    int line_size = cache->line_size;

    // T 2%
    int jlc=0,inca=1,incb=1;
    LineCache *lc;
    while (1) {
      lc = cache->line_cache_start + jlc++;
      lc->linea.resize(line_size);
      lc->lineb.resize(line_size);
      lc->target = location(findx);

      cindx= Indx(nc,1);
      for (int j=0; j<nfree; j++)
	tot_indx[ifree[j]-1] = findx[j];

      int jj=0;
      lc->linear=0;
      while(1) {
	for (int j=0; j<nc; j++) {
	  int k1=icontr[2*j];
	  int k2=icontr[2*j+1];
	  tot_indx[k1-1] = cindx[j];
	  tot_indx[k2-1] = cindx[j];
	}

	// Extract the A and B parts of the indices
	copy(&tot_indx[0],&tot_indx[niA],&aindx[0]);
	copy(&tot_indx[niA],&tot_indx[ndims],&bindx[0]);

	lc->linea[jj] = A.location(aindx);
	lc->lineb[jj] = B.location(bindx);

	if (jj==1) {
	  inca = lc->linea[1] - lc->linea[0];
	  incb = lc->lineb[1] - lc->lineb[0];
	  lc->linear=1;
	} else if (lc->linear && jj>1) {

	  if (lc->linea[jj] - lc->linea[jj-1] != inca) lc->linear=0;
	  if (lc->lineb[jj] - lc->lineb[jj-1] != incb) lc->linear=0;
	}

	jj++;
	if (!inc(cindx,ndimsc)) break;
      }
      lc->starta = &*lc->linea.begin();
      lc->startb = &*lc->lineb.begin();
      lc->inca = inca;
      lc->incb = incb;
      // lc->linear = 0; // force non-linear

      if (!inc(findx,ndimsf)) break;
    }
    // T 97%
    fastmat_stats.tpart += MPI_Wtime()-start;

    // Hay que contar mejor cuantos elementos hay que traer
    // ctx->op_count.get += cache->nlines*cache->line_size;
    int ntot = cache->nlines*cache->line_size;
    ctx->op_count.put += cache->nlines*cache->line_size;
    ctx->op_count.sum += ntot;
    ctx->op_count.mult += ntot;

    psc = new prod_subcache_t(cache);
    assert(psc);
    assert(!cache->sc);
    cache->sc = psc;
    psc->superlinear = 0;
    if (FASTMAT2_USE_DGEMM) {
      psc->ident();
      // if (!psc->not_superlinear_ok()) 
      //   printf("NOT SL!!\n");
      // psc->print();
      
      glob_fm2stats.was_sl_count += psc->superlinear;
      glob_fm2stats.was_not_sl_count += !psc->superlinear;
      
      if (psc->superlinear) {
        glob_fm2stats.was_sl= 1;
        psc->superlinear = glob_fm2stats.use_dgemm;
      }
    }
  }
  fastmat_stats.tcall_not_cached += MPI_Wtime()-start;

  psc = dynamic_cast<prod_subcache_t *> (cache->sc);
  assert(psc);

  if (psc->superlinear) psc->dgemm();
  else {
    // Perform computations using cached addresses
    int 
      nlines = cache->nlines,
      mm=cache->line_size, inca, incb;
    double **pa,**pb,**pa_end,sum,*paa,*pbb,*paa_end;
    LineCache *lc=NULL;

    for (int j=0; j<nlines; j++) {
      lc = cache->line_cache_start+j;
      pa = lc->starta;
      pb = lc->startb;
      inca = lc->inca;
      incb = lc->incb;
      if (lc->linear) {
#if 1
        sum=0.;
        paa = *pa;
        pbb = *pb;
        //         if (inca==1 && incb==1) {
        //           for (int k=0; k<mm; k++) {
        //             sum += (*paa)*(*pbb);
        //             paa++; pbb++;
        //           }
        //         } else 
        if (inca==1) {
          paa_end = paa + mm;
          while (paa<paa_end) {
            sum += (*paa)*(*pbb);
            paa++;
            pbb += incb;
          }
        } else if (incb==1) {
          paa_end = paa + inca*mm;
          while (paa<paa_end) {
            sum += (*paa)*(*pbb);
            paa += inca;
            pbb++;
          }
        } else {
          paa_end = paa + inca*mm;
          while (paa<paa_end) {
            sum += (*paa)*(*pbb);
            paa += inca;
            pbb += incb;
          }
        }
#else
        paa = *pa;
        pbb = *pb;
        paa_end = paa + lc->inca * cache->line_size;
        sum=0.;
        while (paa<paa_end) {
          sum += (*paa)*(*pbb);
          paa += lc->inca;
          pbb += lc->incb;
        }
#endif
      } else {
        pa_end = pa + cache->line_size;
        sum=0.;
        while (pa<pa_end) {
          sum += (**pa++)*(**pb++);
        }
      }
      *(lc->target) = sum;
    }
  }

  if (!ctx->use_cache) delete cache;
  fastmat_stats.tcall += MPI_Wtime()-start;
  return *this;
}
#endif

#ifdef USE_MPROD_FOR_2MATS
#if 0
// OLD VERSION (CALLS PROD)
FastMat2 & 
FastMat2::prod(const FastMat2 & A0,
               const FastMat2 & A1,
               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
    ctx->check(&A0);
    ctx->check(&A1);

    Indx indx;
    indx.push_back(m);
    // READ_INT_ARG_LIST(indx);
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2);
    ctx->check(indx);
  }
#endif

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
    mat_list.push_back(&A0);
    mat_list.push_back(&A1);
    
    Indx &indx = mpwrpsc->indx;
    indx.push_back(m);
    READ_INT_ARG_LIST(indx);
  }

  mpwrpsc = dynamic_cast<mprodwrp_subcache_t *> (cache->sc);
  assert(mpwrpsc);

  prod(mpwrpsc->mat_list,mpwrpsc->indx);
  if (!ctx->use_cache) delete cache;

  return *this;
}
#else
// NEW VERSION (CREATES EXPLICITLY THE MPROD CACHE)
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(const FastMat2 & A0,
               const FastMat2 & A1,
               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
    ctx->check(&A0);
    ctx->check(&A1);

    Indx indx;
    indx.push_back(m);
    // READ_INT_ARG_LIST(indx);
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2);
    ctx->check(indx);
  }
#endif

  FastMatCache *cache = ctx->step();
  multiprod_subcache_t *mpsc=NULL;
  if (!ctx->was_cached) {
    
    assert(!cache->sc);
    vector<const FastMat2 *> mat_list;
    mat_list.push_back(&A0);
    mat_list.push_back(&A1);
    
    Indx indx;
    indx.push_back(m);
    READ_INT_ARG_LIST(indx);
    mpsc = new multiprod_subcache_t(cache,ctx,*this,mat_list,indx);
    cache->sc = mpsc;
  }

  mpsc = dynamic_cast<multiprod_subcache_t *> (cache->sc);
  assert(mpsc);
  mpsc->make_prod();

  if (!ctx->use_cache) delete cache;

  return *this;
}
#endif
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void fastmat_stats_t::print() {
  int ncall_cached = ncall - ncall_not_cached;
  double 
    tcall_cached = tcall - tcall_not_cached,
    nhalf = (tcall_not_cached*ncall_cached)/
    (tcall_cached*ncall_not_cached);
  printf("total calls %d (%g secs, %g rate secs/call)\n"
         "not cached %d, (%g secs, %g rate secs/call)\n"
         "tpart %g (%.2g%%)\n"
         "cached %d, (%g secs, %g rate secs/call)\n"
         "nhalf %f\n",
         ncall,tcall,tcall/ncall,
         ncall_not_cached,tcall_not_cached,tcall_not_cached/ncall_not_cached,
         tpart,tpart/tcall_not_cached*100.0,
         ncall,tcall_cached,tcall_cached/ncall,
         nhalf);
}
