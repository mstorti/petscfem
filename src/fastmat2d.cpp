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
intmax_t 
compute_optimal_order(const mat_info_cont_t &mat_info_cont,
                      vector<int> &order) {

  int qmin,rmin,qfree,rfree,qr1;
  vector<int> ordermin;
  int nmat = mat_info_cont.size();
  intmax_t nopsmin = -1;

  for (int q=0; q<nmat-1; q++) {
    const mat_info &qmi = mat_info_cont[q];
    for (int r=q+1; r<nmat; r++) {
      const mat_info &rmi = mat_info_cont[r];
      // Computes the number of operations
      mat_info smi;
      intmax_t nops = compute_opcount(qmi,rmi,smi,qfree,rfree,qr1);
      vector<int> order1;
      if (nmat>2) {
        mat_info_cont_t mat_info_cont_cpy = mat_info_cont;
        mat_info_cont_cpy[q] = smi;
        mat_info_cont_cpy.erase(mat_info_cont_cpy.begin()+r);
        nops += compute_optimal_order(mat_info_cont_cpy,order1);
      }
      if (nopsmin<0 || nops<nopsmin) {
        qmin = q;
        rmin = r;
        ordermin = order1;
        nopsmin = nops;
      }
    }
  }
  order.clear();
  order.push_back(qmin);
  order.push_back(rmin);

  int n = ordermin.size();
  assert(n==2*(nmat-2));
  for (int j=0; j<n; j++) {
    int oj = ordermin[j];
    assert(oj>=0 && oj<nmat-1);
    order.push_back(oj);
  }

  //#define VERBOSE
#ifdef VERBOSE
  printf("---------------\n");
  printf("nmat %d, nops %jd\n",nmat,nopsmin);
  for (int j=0; j<nmat; j++) {
    const mat_info &mi = mat_info_cont[j];
    printf("a%d: (ctr ",j);
    int rank = mi.contract.size();
    for (int k=0; k<rank; k++)
      printf("%d ",mi.contract[k]);
    printf(") (dims ");
    for (int k=0; k<rank; k++)
      printf("%d ",mi.dims[k]);
    printf(")\n");
  }

  assert(int(order.size())==2*(nmat-1));
#if 0
  for (int j=0; j<nmat-1; j++) 
    printf("prod %d a%d*a%d\n",j,order[2*j],order[2*j+1]);
#endif
#endif

  return nopsmin;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
intmax_t 
compute_heuristic_order(const mat_info_cont_t &mat_info_cont,
                        vector<int> &order) { 

  int nmat = mat_info_cont.size();
  mat_info_cont_t 
    mat_info_cont_cpy = mat_info_cont;
  assert(order.empty());
  int qfree,rfree,qr1;
  mat_info smimin;
  intmax_t nopscount=0;
  for (int j=0; j<nmat-1; j++) {
    int q, r;
    int nmat1=mat_info_cont_cpy.size();
    assert(nmat1==nmat-j);
    int nopsmin=-1;
    for (int jq=0; jq<nmat1-1; jq++) {
      for (int jr=jq+1; jr<nmat1; jr++) {
        mat_info smi;
        const mat_info &qmi = mat_info_cont_cpy[jq];
        const mat_info &rmi = mat_info_cont_cpy[jr];
        int nops = compute_opcount(qmi,rmi,smi,qfree,rfree,qr1);;
        if (nopsmin<0 || nops<nopsmin) {
          nopsmin = nops;
          q = jq;
          r = jr;
          smimin = smi;
        }
      }
    }
    nopscount += nopsmin;
    mat_info_cont_cpy[q] = smimin;
    mat_info_cont_cpy.erase(mat_info_cont_cpy.begin()+r);
    order.push_back(q);
    order.push_back(r);
  }
  assert(int(order.size())==2*(nmat-1));
  return nopscount;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
intmax_t 
compute_natural_order(const mat_info_cont_t &mat_info_cont,
                        vector<int> &order) { 

  int nmat = mat_info_cont.size();
  intmax_t nopscount=0;
  mat_info_cont_t 
    mat_info_cont_cpy = mat_info_cont;
  assert(order.empty());
  int qfree,rfree,qr1,q,r;
  for (int j=0; j<nmat-1; j++) {
#if 1
    // Take the first two
    q=0; r=1;
#else
    // Take the last two 
    r=nmat-j-1;
    q=r-1; 
#endif
    const mat_info &qmi = mat_info_cont_cpy[q];
    const mat_info &rmi = mat_info_cont_cpy[r];
    mat_info smi;
    nopscount += compute_opcount(qmi,rmi,smi,qfree,rfree,qr1);
    mat_info_cont_cpy[q] = smi;
    mat_info_cont_cpy.erase(mat_info_cont_cpy.begin()+r);
    order.push_back(q);
    order.push_back(r);
  }
  assert(int(order.size())==2*(nmat-1));

  return nopscount;
}
