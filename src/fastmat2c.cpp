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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
struct mat_info {
  // Pointers to old matrices should be `const'
  FastMat2 *Ap;
  vector<int> contract;
  vector<int> dims;
  int type;
  mat_info() : Ap(NULL), type(UNKNOWN) {}
  
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static int vfind(vector<int> &v,int x) {
  int n=v.size();
  for (int j=0; j<n; j++)
    if (v[j]==x) return 1;
  return 0;
}

typedef map<int,mat_info> mat_info_cont_t;

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

static int compute_opcount(mat_info &qmi,mat_info &rmi) {
  vector<int> 
    &qc = qmi.contract,
    &qd = qmi.dims,
    &rc = rmi.contract,
    &rd = rmi.dims;
  int 
    rfree=1,
    qfree=1,
    qr1=1,
    qr2=1,
    qrank = qd.size(),
    rrank = rd.size(),
    k,dim;
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// General case
FastMat2 & 
FastMat2::prod(vector<const FastMat2 *> &mat_list,
               Indx &indx) {
  
  // Check sum of ranks of matrices equals
  // size of index vector passed
  int nmat = mat_list.size(), nindx = 0;
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
  mat_info_cont_t mat_info_cont;
  int j=0;
  for (int k=0; k<nmat; k++) {
    FastMat2 &A = *(FastMat2 *)mat_list[k];
    mat_info_cont[k] = mat_info();
    mat_info &m = mat_info_cont[k];
    m.Ap = &A;
    m.type = OLD;
    int rank = A.n();
    m.contract.resize(rank);
    m.dims.resize(rank);
    for (int l=0; l<rank; l++) {
      m.contract[l] = indx[j++];
      m.dims[l] = A.dim(l+1);
    }
  }

  vector<FastMat2 *> new_matrices;
  int new_mat_indx = nmat;
  mat_info_cont_t::iterator q,r,qmin,rmin;
  int qfree,rfree,qr1,qr2,nopsmin,nops,
    qrank,rrank,k,dim,qkey,rkey;
  while (1) {
#if 1
    for (int j=0; j<new_mat_indx; j++) {
      if (mat_info_cont.find(j)!=mat_info_cont.end()) {
        printf("[a%d:",j);
        mat_info &qmi = mat_info_cont.find(j)->second;
        vector<int> &qc = qmi.contract;
        int n=qc.size();
        for (int l=0; l<n-1; l++) printf("%d,",qc[l]);
        if (n>0) printf("%d",qc[n-1]);
        printf("]");
      }
    }
    printf("\n");
#endif
    if (mat_info_cont.size()<=2) break;

    printf("nbr of matrices %zu\n",mat_info_cont.size());
    // search for the product with lowest number
    // of operations
    q = mat_info_cont.begin();
    nopsmin=-1;
    qmin = mat_info_cont.end();
    rmin = mat_info_cont.end();
    while (q != mat_info_cont.end()) {
      qkey = q->first;
      mat_info &qmi = q->second;
      r = q; r++;
      while (r != mat_info_cont.end()) {
        rkey = q->first;
        mat_info &rmi = r->second;
#if 0
        vector<int> 
          &rc = rmi.contract,
          &rd = rmi.dims;
        rfree=1;
        qfree=1;
        qr1=1;
        qr2=1;
        qrank = qd.size();
        rrank = rd.size();
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
#endif
#if 0
        print_mat_info(q,"q: ");
        print_mat_info(r,"r: ");
        
        printf("qfree %d, rfree %d, qr1 %d, qr2 %d \n",
               qfree,rfree,qr1,qr2);
#endif        
        nops = compute_opcount(qmi,rmi);
        if (nopsmin<0 || nops<nopsmin) {
          nopsmin = nops;
          qmin = q;
          rmin = r;
        }
        r++;
      }
      q++;
    }
    assert(nopsmin>0);

    // Insert new entry in 
    int skey = new_mat_indx++;
    mat_info &smi = mat_info_cont[skey];
    smi.type = NEW;
    vector<int> 
      &sc = smi.contract,
      &sd = smi.dims;

    qkey = qmin->first;
    mat_info &qmi = qmin->second;
    vector<int> 
      &qc = qmi.contract,
      &qd = qmi.dims;

    rkey = rmin->first;
    mat_info &rmi = rmin->second;
    vector<int> 
      &rc = rmi.contract,
      &rd = rmi.dims;

    qrank = qd.size();
    rrank = rd.size();
    for (int j=0; j<qrank; j++) {
      k = qc[j];
      dim = qd[j];
      if (k>0 || !vfind(rc,k)) {
        sc.push_back(k);
        sd.push_back(dim);
      }
    }
    for (int j=0; j<rrank; j++) {
      k = rc[j];
      dim = rd[j];
      if (k>0 || !vfind(qc,k)) {
        sc.push_back(k);
        sd.push_back(dim);
      }
    }

#if 0
    mat_info_cont_t::iterator 
      s = mat_info_cont.find(skey);
    print_mat_info(s,"s: ");
#endif

    printf("contracts a%d,a%d\n",qkey,rkey);
    mat_info_cont.erase(qkey);
    mat_info_cont.erase(rkey);
  }

  exit(0);
  return *this;
}
