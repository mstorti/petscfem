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

// //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// // Four matrices involved... The doc should refer to 3 matrices
// FastMat2 & 
// FastMat2::prod(const FastMat2 & A,
//                const FastMat2 & B,
//                const FastMat2 & C,
//                const FastMat2 & D,
//                const int m,INT_VAR_ARGS_ND) {
//   assert(0);
//   return *this;
// }

#include "./mproddef.h"

#define OLD 0
#define NEW 1
#define UNKNOWN 2

#define INACTIVE 0
#define ACTIVE 1
#define UNDEF -1

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
struct mat_info {
  // Pointers to old matrices should be `const'
  FastMat2 *Ap;
  vector<int> contract;
  vector<int> dims;
  int type, is_active, position;
  mat_info() : Ap(NULL), 
               type(UNKNOWN),
               is_active(UNDEF) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static int vfind(vector<int> &v,int x) {
  int n=v.size();
  for (int j=0; j<n; j++)
    if (v[j]==x) return 1;
  return 0;
}

#if 0
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
static int compute_opcount(mat_info &qmi,mat_info &rmi,
                           int &qfree,int &rfree,int &qr1) {
  vector<int> 
    &qc = qmi.contract,
    &qd = qmi.dims,
    &rc = rmi.contract,
    &rd = rmi.dims;
  int 
    qr2=1,
    qrank = qd.size(),
    rrank = rd.size(),
    k,dim;
  rfree=1;
  qfree=1;
  qr1=1;
  for (int j=0; j<qrank; j++) {
    k = qc[j];
    dim = qd[j];
    if (k>0 || !vfind(rc,k)) qfree *= dim;
    else qr1 *= dim;
  }
  for (int j=0; j<rrank; j++) {
    k = rc[j];
    dim = rd[j];
    if (k>0 || !vfind(qc,k)) rfree *= dim;
    else qr2 *= dim;
  }
  assert(qr1==qr2);
  return qfree*rfree*qr1;
}

// Just for debugging makes the first
// contraction that does a nontrivial work
int fastmat_multiprod_use_first=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class multiprod_subcache_t : public FastMatSubCache {
public:
  int nmat;
  vector<mat_info> mat_info_cont;
  vector<int> order;
  multiprod_subcache_t(FastMatCache *cache_a) { }
  ~multiprod_subcache_t();
  void make_prod();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
multiprod_subcache_t::~multiprod_subcache_t() {
  int n = mat_info_cont.size();
  for (int j=0; j<n; j++) {
    mat_info &mi = mat_info_cont[j];
    if (mi.type == NEW && mi.Ap) delete mi.Ap;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void multiprod_subcache_t::make_prod() {
  int nprod = nmat-1,
    qkey, rkey, skey;
  char line[100];
  for (int j=0; j<nprod; j++) {
    mat_info 
      &smi = mat_info_cont[order[3*j]],
      &qmi = mat_info_cont[order[3*j+1]],
      &rmi = mat_info_cont[order[3*j+2]];
    smi.Ap->prod(*qmi.Ap,*rmi.Ap,qmi.contract,rmi.contract);

    sprintf(line,"q[%d]: ",order[3*j+1]);
    qmi.Ap->print(line);

    sprintf(line,"r[%d]: ",order[3*j+2]);
    rmi.Ap->print(line);

    sprintf(line,"s[%d]: ",order[3*j]);
    smi.Ap->print(line);
  }
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
struct compare_position_t : 
  public binary_function<double, double, bool> {
  static vector<mat_info> *mat_info_cont_p;
  bool operator()(double j, double k) { 
    return (*mat_info_cont_p)[j].position 
      < (*mat_info_cont_p)[k].position;
  }
};

vector<mat_info> *compare_position_t::mat_info_cont_p = NULL;
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class active_mat_info_t {
public:
  int key,position;
  bool operator<(const active_mat_info_t y) const {
    if (position != y.position) 
      return position < y.position;
    else 
      return key < y.key;
  }
  active_mat_info_t(int k,int p) 
  : key(k), position(p) { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// General case
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
    mpsc = new multiprod_subcache_t(cache);
    assert(!cache->sc);
    cache->sc = mpsc;
    
    // Check sum of ranks of matrices equals
    // size of index vector passed
    int nmat = mat_list.size(), nindx = 0;
    mpsc->nmat = nmat;
    PETSCFEM_ASSERT0(nmat>=2,"prod needs at least 2 matrices");  

    for (int j=0; j<nmat; j++)
      nindx += mat_list[j]->n();
    PETSCFEM_ASSERT(nindx==indx.size(),"Not correct nbr of indices. "
                    "Sum of ranks %d. nbr of indices %d",nindx,indx.size());  
    // For contracted indices (negeative) count how many times
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
    vector<mat_info> &mat_info_cont = mpsc->mat_info_cont;
    mat_info_cont.resize(nmat);
    vector<int> &order = mpsc->order;
     
    typedef std::set<active_mat_info_t> active_mat_info_cont_t;
    active_mat_info_cont_t active_mat_info_cont;
    int j=0;
    for (int k=0; k<nmat; k++) {
      FastMat2 &A = *(FastMat2 *)mat_list[k];
      mat_info &m = mat_info_cont[k];
      active_mat_info_cont
        .insert(active_mat_info_t(k,k));
      m.Ap = &A;
      m.type = OLD;
      m.position = k;
      int rank = A.n();
      m.contract.resize(rank);
      m.dims.resize(rank);
      for (int l=0; l<rank; l++) {
        m.contract[l] = indx[j++];
        m.dims[l] = A.dim(l+1);
      }
    }

    int new_mat_indx = nmat;
    active_mat_info_cont_t
      ::iterator q,r,qmin,rmin;
    int nact;
    int qfree,rfree,qr1,qr2,nopsmin,nops,
      qrank,rrank,k,dim,qkey,rkey,nopscount=0;
    printf("fastmat_multiprod_use_first %d\n",
           fastmat_multiprod_use_first);
    while (1) {
      nact = active_mat_info_cont.size();
#if 1
      // Print current matrices 
      q = active_mat_info_cont.begin();
      while (q!=active_mat_info_cont.end()) {
        int j = q->key;
        printf("[a%d:",j);
        mat_info &qmi = mat_info_cont[j];
        qmi.is_active = ACTIVE;
        vector<int> &qc = qmi.contract;
        int n=qc.size();
        for (int l=0; l<n-1; l++) printf("%d,",qc[l]);
        if (n>0) printf("%d",qc[n-1]);
        printf("]");
        q++;
      }
      printf("\n");
#endif
      if (active_mat_info_cont.size()<2) break;

      printf("nbr of matrices %zu\n",active_mat_info_cont.size());

      // search for the product with lowest number
      // of operations
      q = active_mat_info_cont.begin();
      nopsmin=-1;
      while (q != active_mat_info_cont.end()) {
        qkey = q->key;
        mat_info &qmi = mat_info_cont[qkey];
        r = q; r++;
        while (r != active_mat_info_cont.end()) {
          rkey = r->key;
          mat_info &rmi = mat_info_cont[rkey];
#if 0
          print_mat_info(q,"q: ");
          print_mat_info(r,"r: ");

          printf("qfree %d, rfree %d, qr1 %d, qr2 %d \n",
                 qfree,rfree,qr1,qr2);
#endif        
          nops = compute_opcount(qmi,rmi,qfree,rfree,qr1);
          if (fastmat_multiprod_use_first?
              qr1>1 
              : (nopsmin<0 || nops<nopsmin)) {
            nopsmin = nops;
            qmin = q;
            rmin = r;
            if (fastmat_multiprod_use_first) break;
          }
          r++;
        }
        q++;
        if (fastmat_multiprod_use_first && nopsmin>0) break;
      }
      assert(nopsmin>0);
      nopscount += nopsmin;

      mat_info_cont.push_back(mat_info());
      // Insert new entry in 
      int skey = new_mat_indx++;
      mat_info &smi = mat_info_cont[skey];
      int reo=1;
      if (nact>2) {
        smi.Ap = new FastMat2;
        smi.type = NEW;
      } else {
        smi.Ap = this;
        smi.type = OLD;
        reo = 0;
      }
      smi.is_active = ACTIVE;
      vector<int> 
        &sc = smi.contract,
        &sd = smi.dims;

      assert(qmin->position<=rmin->position);

      mat_info &qmi = mat_info_cont[qmin->key];
      qmi.is_active = INACTIVE;
      vector<int> 
        &qc = qmi.contract,
        &qd = qmi.dims;

      mat_info &rmi = mat_info_cont[rmin->key];
      rmi.is_active = INACTIVE;
      vector<int> 
        &rc = rmi.contract,
        &rd = rmi.dims;

      qrank = qd.size();
      rrank = rd.size();
      int sindx = 1;
      for (int j=0; j<qrank; j++) {
        k = qc[j];
        dim = qd[j];
        if (k>0 || !vfind(rc,k)) {
          if (reo) qc[j] = sindx++;
          sc.push_back(k);
          sd.push_back(dim);
        }
      }
      for (int j=0; j<rrank; j++) {
        k = rc[j];
        dim = rd[j];
        if (k>0 || !vfind(qc,k)) {
          if (reo) rc[j] = sindx++;
          sc.push_back(k);
          sd.push_back(dim);
        }
      }
      order.push_back(skey);
      smi.position = qmi.position;
      order.push_back(qmin->key);
      order.push_back(rmin->key);

      active_mat_info_cont.erase(qmin);
      active_mat_info_cont.erase(rmin);
      active_mat_info_cont
        .insert(active_mat_info_t(skey,smi.position));

      int n=order.size();
      printf("contracts a%d,a%d\n",qmin->key,rmin->key);
    }
    printf("total ops count %d\n",nopscount);
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

      map<int,int> new_neg;
      int new_neg_free=-1;
      int k,
        qrank = qc.size(),
        rrank = rc.size();
      for (int l=0; l<qrank; l++) {
        k = qc[l];
        if (k<0) {
          if (new_neg.find(k)==new_neg.end()) 
            new_neg[k] = new_neg_free--;
          qc[l] = new_neg[k];
        }
      }
      for (int l=0; l<rrank; l++) {
        k = rc[l];
        if (k<0) {
          if (new_neg.find(k)==new_neg.end()) 
            new_neg[k] = new_neg_free--;
          rc[l] = new_neg[k];
        }
      }
#if 1
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

  mpsc = dynamic_cast<multiprod_subcache_t *> (cache->sc);
  assert(mpsc);
  mpsc->make_prod();

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
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

  prod_subcache_t *psc=NULL;
  if (!ctx->was_cached  ) {
    Indx ia,ib,ii;

    Indx Afdims,Bfdims,fdims;
    A.get_dims(Afdims);
    B.get_dims(Bfdims);

    int niA = Afdims.size();
    int niB = Bfdims.size();
    int ndims = niA+niB;
    
    for (unsigned int j=0; j<ixa.size(); j++) 
      ii.push_back(ixa[j]);
    for (unsigned int j=0; j<ixb.size(); j++) 
      ii.push_back(ixb[j]);

    int nfree=0,nc=0;
    for (int j=0; j<ii.size(); j++) {
      int k = ii[j];
      if (k>0 && k>nfree) nfree = k;
      if (k<0 && -k>nc) nc = -k;
    }

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
    if (comp_storage_size(ndimsf)==0) {
      cache->nelems=0;
      return *this;
    }

    get_dims(fdims);
    if (ndimsf != fdims) {
      Afdims.print("A free dims: ");
      Bfdims.print("B free dims: ");
      ndimsf.print("Combined free dims: ");
      fdims.print("Free dims on result matrix: ");
      printf("Combined free dims doesn't match free"
	     " dims of  result.\n");
      abort();
    }
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
      assert(nd1==nd2);
      ndimsc[j]=nd1;
    }

    // ndimsc.print("dimensions of the contracted part: ");

    Indx findx(nfree,1),cindx(nc,1),tot_indx(ndims,0),
      aindx(niA,0),bindx(niB,0);

  // Loading addresses in cache
  // For each element in the distination target, we store the complete
  // list of addresses of the lines of elements that contribute to
  // it. 
    cache->nlines = mem_size(ndimsf);
    cache->prod_cache.resize(cache->nlines);
    cache->line_cache_start = &*cache->prod_cache.begin();
    cache->line_size = mem_size(ndimsc);
    int line_size = cache->line_size;

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
    psc->ident();
    if (FASTMAT2_USE_DGEMM) {
      if (!psc->not_superlinear_ok()) {
        printf("NOT SL!!\n");
      }
      // psc->print();
      
      glob_fm2stats.was_sl_count += psc->superlinear;
      glob_fm2stats.was_not_sl_count += !psc->superlinear;
      
      if (psc->superlinear) {
        glob_fm2stats.was_sl= 1;
        psc->superlinear = glob_fm2stats.use_dgemm;
      }
    }
  }

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
  return *this;
}

#ifdef USE_MPROD_FOR_2MATS
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(const FastMat2 & A0,
               const FastMat2 & A1,
               const int m,INT_VAR_ARGS_ND) {

  vector<const FastMat2 *> mat_list;
  mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  
  Indx indx;
  indx.push_back(m);
  READ_INT_ARG_LIST(indx);
  prod(mat_list,indx);
  return *this;
}
#endif
