// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: project.h,v 1.1 2005/03/02 23:26:31 mstorti Exp $
#ifndef PETSCFEM_PROJECT_H
#define PETSCFEM_PROJECT_H

class FemInterp {
private:
  dvector<double> xnod;
  dvector<int> icone;

  ANNkd_tree *kdtree;
  vector<ANNidx> nn_idx_v;
  vector<ANNdist> nn_dist_v;
  // ANNpoint nn;
  ANNpointArray pts;
  
  int knbr,ndim,nnod,ndimel,
    nel,nelem,ndof,nd1;

  FastMat2 C,C2,
    invC,invC2, invCt,
    x2,dx2, x2prj,x2prjmin, x12,
    x13,x1, nor,L, b,u1_loc, u2,
    Lmin;

  FastMatCachePosition cp,cp1;
  FastMatCacheList cache_list;

  vector<int> restricted;
public:
  int use_cache;
  double tol;

  FemInterp();

  ~FemInterp();

  void clear();

  // `knbr_a' = nbr of neighbors to take (set to 10 usually)
  // `ndof_a' nbr of degrees of reedom
  // `ndimel_a' dimension of the manifold
  // 
  void init(int knbr_a, int ndof_a, int ndimel_a,
	    const dvector<double> &xnod_a,
	    const dvector<int> &icone_a);

  void interp(const dvector<double> &xnod2,
	      const dvector<double> &u,
	      dvector<double> &ui);
};


#endif
