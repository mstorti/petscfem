//__INSERT_LICENSE__
#include <cmath>
#include <ctime>
#include <cstdio>
#include <stdint.h>
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
intmax_t compute_optimal_order(const mat_info_cont_t &mat_info_cont,
                               vector<int> &order) {

  int qmin,rmin,qfree,rfree,qr1;
  vector<int> ordermin;
  int nmat = mat_info_cont.size();
  intmax_t nopsmin = -1;
  

  for (int q=0; q<nmat; q++) {
    const mat_info &qmi = mat_info_cont[q];
    for (int r=q+1; q<nmat; q++) {
      const mat_info &rmi = mat_info_cont[r];
      // Computes the number of operations
      mat_info smi;
      intmax_t nops = compute_opcount(qmi,rmi,smi,qfree,rfree,qr1);
      vector<int> order1;
      if (nmat>2) {
        mat_info_cont_t mat_info_cont_cpy = mat_info_cont;
        mat_info_cont_cpy.erase(mat_info_cont_cpy.begin()+r);
        mat_info_cont_cpy.erase(mat_info_cont_cpy.begin()+q);
        mat_info_cont_cpy.push_back(smi);
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

  for (int j=0; j<n; j++) 
    order.push_back(ordermin[j]);
  return nopsmin;
}
