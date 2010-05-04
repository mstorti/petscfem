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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & FastMat2::prod(const FastMat2 & A,
                          const FastMat2 & B,
                          const FastMat2 & C,
                          const int m,INT_VAR_ARGS_ND) {

  vector<const FastMat2 *> mat_list;
  mat_list.push_back(&A);
  mat_list.push_back(&B);
  mat_list.push_back(&C);
  
  Indx indx;
  indx.push_back(m);
  READ_INT_ARG_LIST(indx);
  prod(mat_list,indx);
  return *this;
}

// Four matrices involved... The doc should refer to 3 matrices
FastMat2 & 
FastMat2::prod(const FastMat2 & A,
               const FastMat2 & B,
               const FastMat2 & C,
               const FastMat2 & D,
               const int m,INT_VAR_ARGS_ND) {
  assert(0);
  return *this;
}

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

  map<int,int>::iterator q = indices.begin();
  int erro=0, outrank=-1;
  while (q!=indices.end()) {
    int k=q->first, n=q->second;
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
    q++;
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

  vector<vector<int> > mat_indx(nmat);
  int j=0;
  for (int k=0; k<nmat; k++) {
    int rank = mat_list[k]->n();
    vector<int> &v = mat_indx[k];
    for (int l=0; l<rank; l++) 
      v.push_back(indx[j++]);
  }
  exit(0);
  return *this;
}
