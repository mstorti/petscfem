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
void compute_optimal_order(active_mat_info_cont_t 
                           &active_mat_info_cont) {
#if 0
  active_mat_info_cont_t
    ::iterator q,r,qmin,rmin;
  // nbr of active matrices
  int nact;
  int qfree,rfree,qr1,qr2,nopsmin,nops,
    qrank,rrank,k,dim,qkey,rkey,nopscount=0;
#ifdef VERBOSE
    printf("fastmat_multiprod_use_first %d\n",
           fastmat_multiprod_use_first);
#endif
    while (1) {
      // For each execution of this loop we
      // look for the lowest opcount. The the
      // two corresponding matrices become inactive
      // and the new temporary is created. 

      // Nbr of current active matrices
      nact = active_mat_info_cont.size();

#ifdef VERBOSE
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

      // printf("nbr of matrices %zu\n",active_mat_info_cont.size());

      // search for the product with lowest number
      // of operations
      // q,r, iterate, qmin,rmin keep the lowest combination
      q = active_mat_info_cont.begin();
      nopsmin=-1;
      while (q != active_mat_info_cont.end()) {
        // we identify the matrix by his position
        // in `mat_info_cont' (the `key')
        qkey = q->key;
        mat_info &qmi = mat_info_cont[qkey];
        r = q; r++;
        while (r != active_mat_info_cont.end()) {
          // Same for `r'
          rkey = r->key;
          mat_info &rmi = mat_info_cont[rkey];
#ifdef VERBOSE 
          print_mat_info(q,"q: ");
          print_mat_info(r,"r: ");

          printf("qfree %d, rfree %d, qr1 %d, qr2 %d \n",
                 qfree,rfree,qr1,qr2);
#endif        
          // Computes the number of operations
          nops = compute_opcount(qmi,rmi,qfree,rfree,qr1);
          // `fastmat_multiprod_use_first' is just for debugging.
          // If activated it corresponds to do the first product
          // that does a contraction. For instance in a series
          // of products of 2-rank matrices a=b*c*d*e*...
          // it corresponds to do a=(((b*c)*d)*e)*...
          // Beware it may fail. in case there is no contraction!!
          // For instance a.prod(b,c,1,2);
          if (fastmat_multiprod_use_first?
              (nops>0 || qr1>1)
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
      // Keeps the total number of operations performed
      nopscount += nopsmin;

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
        smi.Ap = new FastMat2(this->ctx);
        smi.type = TMP;
      } else {
        // In this case the output matrix is *this
        smi.Ap = this;
        // It is NOT a temporary
        smi.type = OLD;
        // Do not reorder the indices
        reo = 0;
      }
      // This new matrix is active
      smi.is_active = ACTIVE;
      vector<int> 
        &sc = smi.contract,
        &sd = smi.dims;

      // Check that the active matrices are ordered FIRST
      // by position 
      assert(qmin->position<=rmin->position);

      // Product to be done is S = Q * R
      // get info for Q matrix
      mat_info &qmi = mat_info_cont[qmin->key];
      qmi.is_active = INACTIVE;
      vector<int> 
        &qc = qmi.contract,
        &qd = qmi.dims;

      // get info for R matrix
      mat_info &rmi = mat_info_cont[rmin->key];
      rmi.is_active = INACTIVE;
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
      order.push_back(qmin->key);
      order.push_back(rmin->key);
      
      if (fastmat_stats.print_prod_order) 
        printf("will do a%d = a%d*a%d\n",
               smi.position,
               qmi.position,rmi.position);

      // Mark the Q and R matrices as inactive,
      // and the result S as active
      active_mat_info_cont.erase(qmin);
      active_mat_info_cont.erase(rmin);
      active_mat_info_cont
        .insert(active_mat_info_t(skey,smi.position));

#ifdef VERBOSE
      int n=order.size();
      printf("contracts a%d,a%d\n",qmin->key,rmin->key);
#endif
    }

#endif

}
