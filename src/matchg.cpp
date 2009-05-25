#include <vector>
#include <algorithm>

#include <src/utils.h>
#include <src/dvector.h>
#include <src/dvector2.h>

extern "C" {
#define __log2 ___log2
#define drand48 __drand48
#define srand48 __srand48
#include <metis.h>
#undef __log2
#undef __drand48
#undef __srand48
}

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void make_sym_graph(dvector<double> &G) {
  PETSCFEM_ASSERT0(G.rank()==2,"Graph G must have rank 2");  
  int N = G.size(0);
  PETSCFEM_ASSERT0(G.size(1)==N,"Graph G must be a square matrix");  
  for (int j=0; j<N-1; j++) {
    G.e(j,j) = 0.0;
    for (int k=j+1; j<N; j++) {
      double w = maxd(2,G.e(j,k),G.e(k,j));
      G.e(j,k) = w;
      G.e(k,j) = w;
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int is_sym_graph(const dvector<double> &G, 
                   double tol=1e-10) {
  PETSCFEM_ASSERT0(G.rank()==2,"Graph G must have rank 2");  
  int N = G.size(0);
  PETSCFEM_ASSERT0(G.size(1)==N,"Graph G must be a square matrix");  
  for (int j=0; j<N-1; j++) {
    if (G.e(j,j)!=0.0) return 0;
    for (int k=j+1; j<N; j++) {
      if (fabs(G.e(j,k)-G.e(k,j))>tol) {
        return 0;
      }
    }
  }
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Find a mapping between domains and processors
// so that to pair the data-flux matrix (bflux)
// with bandwidth matrix `bw'.
// `bflux(dj,dk)' is the amount of communication that
// is expected between domains `dj' and `dk'. 
// `bw(pj,pk)' is the bandwidth between processors
// pj and pk. We look for a mapping p = proc(d), or
// the inverse d = domain(p), that minimizes the rate
// max (k>j) bflux(j,k)/bw(j,k). 
void match_graph(const dvector<double> &bw,
                 const dvector<double> &bflux,
                 dvector<int> &proc) {
  //#define DBG
  int N = proc.size();
  PETSCFEM_ASSERT0(proc.rank()==1,"error in proc rank");  

  PETSCFEM_ASSERT0(bw.rank()==2,"error in bw rank");  
  PETSCFEM_ASSERT0(bw.size(0)==N,"error in bw dim 0");  
  PETSCFEM_ASSERT0(bw.size(1)==N,"error in bw dim 1");  

  PETSCFEM_ASSERT0(bflux.rank()==2,"error in bflux rank");  
  PETSCFEM_ASSERT0(bflux.size(0)==N,"error in bflux dim 0");  
  PETSCFEM_ASSERT0(bflux.size(1)==N,"error in bflux dim 1");  

  proc.set(-1);
//   PETSCFEM_ASSERT0(is_sym_graph(bw),
//                    "Bandwidth graph `bw' must be symmetric");  
//   PETSCFEM_ASSERT0(is_sym_graph(bflux),
//                    "Data flux graph `bflux' must be symmetric");  
  
  dvector<int> domain;
  domain.mono(N).reshape(1,N).set(-1);

  double max_flux=NAN;
  int dj=-1, dk=-1;
  // Compute next max flux to be assigned
  for (int j=0; j<N; j++) {
    for (int k=0; k<N; k++) {
      if (j==k || (proc.e(j)>=0 && proc.e(k)>=0)) continue;
      if (isnan(max_flux) || bflux.e(j,k)>max_flux) {
        max_flux = bflux.e(j,k);
        dj = j; dk = k;
      }
    }
  }
#ifdef DBG
  printf("max flux edge (%d,%d) -> %f\n",dj,dk,max_flux);
#endif

  // Find next channel with max bandwidth
  double max_bw = NAN;
  int pj=-1, pk=-1;
  for (int j=0; j<N; j++) {
    for (int k=0; k<N; k++) {
      if (j==k || (proc.e(j)>=0 && proc.e(k)>=0)) continue;
      if (isnan(max_bw) || bw.e(j,k)>max_bw) {
        max_bw = bw.e(j,k);
        pj = j; pk = k;
      }
    }
  }

  // Assign channel pj,pk to dj,dk
  proc.e(dj) = pj;
  proc.e(dk) = pk;
  domain.e(pj) =dj;
  domain.e(pk) =dk;
#ifdef DBG
  printf("assigning dom=%d -> proc=%d\n",dj,pj);
  printf("assigning dom=%d -> proc=%d\n",dk,pk);
#endif

  while(1) {
    // Find next channel with one vertex assigned and
    // max flux
    max_flux = NAN;
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
        if ((proc.e(j)>=0) == (proc.e(k)>=0)) continue;
        if (isnan(max_flux) || bflux.e(j,k)>max_flux) {
          max_flux = bflux.e(j,k);
          dj = j; dk = k;
        }
      }
    }
    if (isnan(max_flux)) break;
#ifdef DBG
    printf("max flux edge (%d,%d) -> %f\n",dj,dk,max_flux);
#endif

    // Which vertex is assigned
    int asj = (proc.e(dj)>=0);
    int ask = (proc.e(dk)>=0);

    // Verify not both are assigned
    assert(asj != ask);
    // We set `df' to the domain that has a processor
    // already assigned and `dnf' to the other. 
    int df = dj, dnf = dk;
    if (ask) { df = dk, dnf = dj; }
    double max_bw = NAN;
    // `pf' is the processor assigned to `df'
    int pf = proc.e(df), pnf=-1;
    // Find `pnf' as the processor with max bw to `pf'
    assert(pf>=0);
    for (int p=0; p<N; p++) {
      if (domain.e(p)>=0) continue;
      if (isnan(max_bw) || bw.e(pf,p)>max_bw) {
        max_bw = bw.e(pf,p);
        pnf = p;
      }
    }
    assert(!isnan(max_bw));
    proc.e(dnf) = pnf;
    domain.e(pnf) = dnf;
#ifdef DBG
    printf("assigning dom=%d -> proc=%d\n",dnf,pnf);
#endif
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void repart(const dvector<double> &bflux,
            dvector<int> &proc, int ncore) {

  int ndom = proc.size();
  dvector<int> xadj,adjncy,adjwgt, dpart;
  int edgew_scale = 100000;
  PETSCFEM_ASSERT_GE0(bflux.rank()==2,"Internal error\n");
  PETSCFEM_ASSERT_GE0(bflux.size(0)==ndom,"Internal error\n");
  PETSCFEM_ASSERT_GE0(bflux.size(1)==ndom,"Internal error\n");

  double bflux_max = NAN;
  for (int j=0; j<ndom; j++) {
    for (int k=0; k<ndom; k++) {
      if (j==k) continue;
      if (isnan(bflux_max) || bflux.e(j,k)>bflux_max)
        bflux_max = bflux.e(j,k);
    }
  }

  int nedges = ndom*(ndom-1);
  xadj.mono(ndom+1).set(0);
  adjncy.mono(nedges).set(0);
  adjwgt.mono(nedges).set(0);
  dpart.mono(ndom).set(0);
  
  int indx = 0;
  for (int j=0; j<ndom; j++) {
    for (int k=0; k<ndom; k++) {
      if (k==j) continue;
      adjncy.e(indx) = k;
      adjwgt.e(indx) = int(edgew_scale*bflux.e(j,k)/bflux_max);
      if (adjwgt.e(indx)==0) adjwgt.e(indx)=1;
      indx++;
    }
    xadj.e(j+1) = indx;
  }

  int wgtflag=1, numflag=0;
  int nparts = ndom/ncore, options=0, edgecut;
  
  METIS_PartGraphKway(&ndom,xadj.buff(),adjncy.buff(),NULL, 
                      adjwgt.buff(),&wgtflag,&numflag,&nparts, 
                      &options,&edgecut,dpart.buff());

  // Number first all domains in dpart=0,
  // then dpart=1, etc...
  int p=0;
  for (int j=0; j<nparts; j++) 
    for (int k=0; k<ndom; k++) 
      if (dpart.e(k)==j) proc.e(k)=p++;
  assert(p==ndom);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Compute the efficiency of the matching,
// perfo_max = maximum { bflux(dj,dk)/bflux(pj,pk) }
// perfo_sum = sum { bflux(dj,dk)/bflux(pj,pk) }
void perfo(const dvector<double> &bw,
           const dvector<double> &bflux,
           const dvector<int> &proc, double &perfo_max,
           double &perfo_sum) {
  perfo_max = NAN;
  perfo_sum = 0.0;
  int N = proc.size();
  for (int dj=0; dj<N; dj++) {
    for (int dk=0; dk<N; dk++) {
      if (dj==dk) continue;
      int pj = proc.e(dj);
      int pk = proc.e(dk);
      double rate = bflux.e(dj,dk)/bw.e(pj,pk);
      if (isnan(perfo_max) || rate>perfo_max) 
        perfo_max = rate;
      perfo_sum += rate;
    }
  }
}
