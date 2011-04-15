//__INSERT_LICENSE__
// $Id$

#include <cstdio>
#include <cmath>
#include <ctime>
#include <stdint.h>
#include <set>
#include <algorithm>

#include <mpi.h>

#include <src/fastmat2.h>
#include <mkl_cblas.h>
#include <src/fm2stats.h>
#include <src/fem.h>
#include <src/dvector.h>
#include <src/dvector2.h>

double drand1(double) {
  return drand();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void test() {
  int Nmax=100,Nmin=3,Nmat,Nmat1=10;
  int all_ok = 1;
  double tol=1e-5;
  for (Nmat=2; Nmat<=Nmat1; Nmat++) {
    FastMat2::CacheCtx2 ctx;
    // ctx.mprod_order = FastMat2::CacheCtx::mixed;
    vector<const FastMat2 *> mats(Nmat);
    vector<int> dims(Nmat+1);
    for (int j=0; j<Nmat+1; j++) {
      double D = exp(log(Nmin)+drand()*(log(Nmax)-log(Nmin)));
      dims[j] = int(round(D));
      // printf("j %d, dim %d\n",j,dims[j]);
    }
    Indx indx;

    for (int j=0; j<Nmat; j++) {
      FastMat2 *ap = new FastMat2(&ctx,2,dims[j],dims[j+1]);
      ap->fun(drand1);
      mats[j] = ap;
      indx.push_back(-j);
      indx.push_back(-(j+1));
    }
    indx[0] = 1;
    indx[2*Nmat-1] = 2;
    
    intmax_t nopscount = 0;
    for (int j=0; j<Nmat-1; j++) 
      nopscount += 2*dims[0]*dims[j+1]*dims[j+2];

    int ntimes=100;
    FastMat2 res(&ctx);
    
    FastMat2::CacheCtx2::Branch br;
    // ctx.mprod_order = FastMat2::CacheCtx2::natural;
    ctx.use_cache = 1;
    double start = MPI_Wtime();
    for (int j=0; j<ntimes; j++) {
      ctx.jump(br);
      res.prod(mats,indx);
    }
    double elapsed = MPI_Wtime()-start;
    double rate = nopscount*ntimes/elapsed*1.0e-9;
    printf("Nmat %d, nopscount %jd, ntimes %d, "
           "elapsed %f secs, rate %f Gflops\n",
           Nmat, nopscount, ntimes, elapsed, rate);
    // res.print("using multiprod: ");
    
    // Compute using standard 2-matrix product
    FastMat2 res2(&ctx),tmp(&ctx);
    ctx.use_cache = 0;
    ctx.was_cached = 0;
    tmp.set(*mats[0]);
    for (int j=0; j<Nmat-1; j++) {
      res2.clear();
      res2.prod(tmp,*mats[j+1],1,-1,-1,2);
      tmp.clear().set(res2);
    }
    // res2.print("using 2-nat prod\n ");
    res2.minus(res);
    double erro = res2.norm_2_all();
    printf("Nmat %d: error %g, ok<tol? %d\n",Nmat,erro,erro<tol);
    all_ok &= erro<tol;
    for (int j=0; j<Nmat; j++) {
      FastMat2 *ap = (FastMat2 *)mats[j];
      delete ap;
    }
  }
    printf("All OK? %d\n",all_ok);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int main(int argc,char**args) {
  PetscFemInitialize(&argc,&args,NULL,NULL);
  test();
  PetscFinalize();
  return 0;
}
