//__INSERT_LICENSE__
// $Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

#include <cstdio>
#include <string>
#include <mpi.h>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <ANN/ANN.h>

#include "./project.h"

string print_area_coords;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FemInterp::clear() {
  if (kdtree) delete kdtree;
  kdtree = NULL;
  nn_dist_v.clear();
  if (pts) annDeallocPts(pts);
  pts = NULL;
  nn_idx_v.clear();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FemInterp::FemInterp() : 
  kdtree(NULL), pts(NULL), 
  C(&ctx),C2(&ctx),invC(&ctx),invC2(&ctx), invCt(&ctx),
  x2(&ctx),dx2(&ctx), x2prj(&ctx),x2prjmin(&ctx), x12(&ctx),
  x13(&ctx),x1(&ctx), nor(&ctx),L(&ctx), b(&ctx),
  u1_loc(&ctx), u2(&ctx), Lmin(&ctx),
  use_cache(1), use_delaunay(0), 
  tol(1e-6) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FemInterp::~FemInterp() { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FemInterp::init(int knbr_a, int ndof_a, int ndimel_a,
		     const dvector<double> &xnod_a,
		     const dvector<int> &icone_a) {
    
  clear();

  ndof = ndof_a;
  knbr = knbr_a;
  ndimel = ndimel_a;

  // nn_idx = new ANNidx[knbr];
  // nn_dist = new ANNdist[knbr];
  // nn = annAllocPt(ndim);
    
  nnod = xnod_a.size(0);
  ndim = xnod_a.size(1);
  nelem = icone_a.size(0);
  nel = icone_a.size(1);
  brv1.init(nel); 
  brv2.init(nel);

  xnod.clone(xnod_a);
  icone.clone(icone_a);

  // Build ANN octree
  printf("computing element centers ...\n");
  double start = MPI_Wtime();
  FastMat2::CacheCtx2 ctx1;
  ctx1.use_cache = 1;
  FastMat2::CacheCtx2::Branch br1;
  FastMat2 xe(&ctx1,1,ndim),xn(&ctx1,1,ndim);
  if (!use_delaunay) {
    double inel = 1./nel;
    pts = annAllocPts(nelem,ndim);
    // fixme:= should `pts' be freed after??
    // seems that no
    for (int k=0; k<nelem; k++) {
      ctx1.jump(br1);
      xe.set(0.);
      for (int j=0; j<nel; j++) {
        int node = icone.e(k,j);
        xn.set(&xnod.e(node-1,0));
        xe.add(xn);
      }
      xe.scale(inel);
#if 0
      printf("elem %d, x ",k);
      xe.print();
#endif
      for (int j=0; j<ndim; j++)
        pts[k][j] = xe.get(j+1);
    }
    printf("end computing element centers... elapsed %f\n",
           MPI_Wtime()-start);

    printf("loading points in ann-tree ...\n");
    start = MPI_Wtime();
    kdtree = new ANNkd_tree(pts,nelem,ndim);
    printf("end loading points in ann-tree... elapsed %f\n",
           MPI_Wtime()-start);
  } else {
    pts = annAllocPts(nnod,ndim);
    for (int k=0; k<nnod; k++) {
      for (int j=0; j<ndim; j++)
        pts[k][j] = xnod.e(k,j);
    }
    printf("loading points in ann-tree ...\n");
    start = MPI_Wtime();
    kdtree = new ANNkd_tree(pts,nnod,ndim);
    printf("end loading points in ann-tree... elapsed %f\n",
           MPI_Wtime()-start);
  }

  nd1 = ndim+1;
  C.resize(2,nd1,nd1);
  C2.resize(2,nd1,nd1);
  invC.resize(2,nd1,nd1);
  invC2.resize(2,nd1,nd1);
  invCt.resize(2,nd1,nd1);
  x2.resize(1,ndim);
  dx2.resize(1,ndim);
  x2prj.resize(1,ndim);
  x2prjmin.resize(1,ndim);
  x12.resize(1,ndim);
  x13.resize(1,ndim);
  x1.resize(1,ndim);
  nor.resize(1,ndim);
  L.resize(1,nd1);
  Lmin.resize(1,nd1);
  b.resize(1,ndim+1);
  u1_loc.resize(2,ndof,nel);
  u2.resize(1,ndof);

  restricted.resize(nel,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FemInterp::init(int knbr, int ndof, int ndimel,
		     const dvector<double> &xnod) {
  use_delaunay = 0;
  dvector<int> icone;
  icone.reshape(2,0,nel);
  init(knbr, ndof, ndimel,xnod,icone);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FemInterp::interp(const dvector<double> &xnod2,
		       const dvector<double> &u,
		       dvector<double> &ui) {

  double tryav=0.0;
  int nnod2 = xnod2.size(0);
  assert(xnod2.size(1) == ndim);
  ui.a_resize(2,nnod2,ndof).defrag();

  // Minimum distance to elements in mesh1
  double d2min=NAN;
  // Element in mesh1 with minimum distance
  int k1min=0; 

  nn_idx_v.resize(knbr);
  ANNidx *nn_idx = &nn_idx_v[0];
  nn_dist_v.resize(knbr);
  ANNdist *nn_dist = &nn_dist_v[0];
  ANNpoint nn = annAllocPt(ndim);
  // If `ndimel==ndim' then normally we find a point
  // a zero distance of the test point. If `ndimel<ndim'
  // then this is rarely the case and we check only for
  // the `knbr' elements reported by `ANN'. 
  // int nelem_check = (ndimel==ndim? nelem+knbr : knbr);
  int nelem_check = knbr;
  // printf("start interpolation...\n");
  double start = MPI_Wtime();
  FILE *fid=NULL;
  if (print_area_coords.size()>0) 
    fid = fopen(print_area_coords.c_str(),"w");
  ctx.use_cache = use_cache;
  for (int n2=0; n2<nnod2; n2++) {
    ctx.jump(br1);
    x2.set(&xnod2.e(n2,0));
    for (int j=0; j<ndim; j++) 
      nn[j] = xnod2.e(n2,j);
    kdtree->annkSearch(nn,knbr,nn_idx,nn_dist,0.0);
#if 0
    for (int k=0; k<knbr; k++)
      printf("(%d %f) ",nn_idx[k],nn_dist[k]);
    printf("\n");
#endif
    int q;

    for (q=0; q<nelem_check; q++) {
      ctx.jump(br2);
      
      int k = (q<knbr ? nn_idx[q] : q-knbr);
      C.is(1,1,ndim);
      for (int j=0; j<nel; j++) {
	int node = icone.e(k,j);
	C.ir(2,j+1).set(&xnod.e(node-1,0));
      }
      C.rs();
      C.ir(1,nd1).is(1,1,nel).set(1.0).rs();
      // If `ndimel<ndim' complete columns to `ndim+1' with 
      // vectors normal to the element surface, and set the
      // elements in the `ndim+1' row to 0. 
      
      b.is(1,1,ndim).set(x2).rs();
      b.setel(1.0,nd1);
      for (int j=0; j<nel; j++)
	restricted[j] = 0;
      invC.inv(C);
      C2.set(C);
      // b.print("b");
      // C.print("C");
      invCt.ctr(invC,2,1).scale(-1.0);
      invCt.ir(1,nd1).set(0.).rs();
      // invCt.print("invCt");
      int iter=0;
      // FastMat2::branch();
      while(true) {
	ctx.jump(br3);
	// FastMat2::choose(iter);
	iter++;
	invC2.inv(C2);
	L.prod(invC2,b,1,-1,-1);
	// C2.print("C2");
	// L.print("L");
	int neg=0;
        double *Lp = L.storage_begin();
	for (int j=1; j<=nel; j++) {
	  double lj = Lp[j-1];
	  if (lj<-tol) {
	    neg=1;
	    int &flag = restricted[j-1];
	    flag = !flag;
	    if (flag) {
              ctx.jump(brv1(j-1));
	      C2.ir(2,j);
	      invCt.ir(2,j);
	      C2.set(invCt).rs();
	      invCt.rs();
	    } else {
              ctx.jump(brv2(j-1));
	      C2.ir(2,j);
	      C.ir(2,j);
	      C2.set(C).rs();
	      C.rs();
	    }
	  }
	}
	if (!neg) break;
	// assert(iter<=ndim);
	int nitmax=20;
	if (iter>nitmax) {
	  printf("failed to converge in %d iters\n",nitmax);
	  break;
	}
      }
      ctx.jump(br7);

      // Set area coordinates for restricted
      // nodes to 0. 
      double *Lp = L.storage_begin();
      for (int j=0; j<nel; j++) 
	if (restricted[j]) Lp[j] = 0;

      for (int j=nel; j<nd1; j++) Lp[j] = 0;

      // Set area coordinates for normal
      // vectors to 0.
      C.is(1,1,ndim);
      x2prj.prod(C,L,1,-1,-1);
      C.rs();
      // Form distance vector
      dx2.set(x2).minus(x2prj);
      double d2 = dx2.norm_p_all(2.0);
      if (q==0 || d2<d2min) {
        ctx.jump(br8);
	d2min = d2;
	k1min = k;
	Lmin.set(L);
	x2prjmin.set(x2prj);
      }
      if (d2min<tol) break;
    }
    ctx.jump(br9);
    // x2.print("x2");
    // printf("tries %d\n",q+1);
    // x2prjmin.print("x2prjmin");
    // Load state values in `u1_loc'
    for (int j=1; j<=nel; j++) {
      u1_loc.ir(2,j);
      int node = icone.e(k1min,j-1);
      u1_loc.set(&u.e(node-1,0));
    }
    if (print_area_coords.size()>0) {
      fprintf(fid,"%d %d ",n2+1,k1min+1);
      const double *L = Lmin.storage_begin();
      for (int j=0; j<nel; j++)
        fprintf(fid,"%.15f ",L[j]);
      fprintf(fid,"\n");
    }
    u1_loc.rs();
    // Interpolate
    Lmin.is(1,1,nel);
    u2.prod(u1_loc,Lmin,1,-1,-1);
    Lmin.rs();
    tryav += q+1;
    u2.export_vals(&ui.e(n2,0));
  }
  if (print_area_coords.size())
    fclose(fid);
  // printf("end interpolation... elapsed %f\n",MPI_Wtime()-start);
  annDeallocPt(nn);
  // delete[] nn_idx;
  // printf("Averg. nbr of tries %f\n",tryav/nnod2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void nod_vol(const dvector<double> &xnod,
	     const dvector<int> &icone,
	     dvector<double> &area) {
  int ndim = xnod.size(1);
  int nel = icone.size(1);
  int ndimel = nel-1; // Assume that the elements are simplices
  assert(ndimel==ndim-1);
  FastMat2 a(1,ndim),b(1,ndim),c(1,ndim),
    x1(1,ndim), x2(1,ndim), x3(1,ndim);
  int nnod = xnod.size(0);
  area.a_resize(1,nnod).defrag().set(0.);
  int nelem = icone.size(0);
  double area_tot = 0.0;
  for (int j=0; j<nelem; j++) {
    int n1 = icone.e(j,0);
    int n2 = icone.e(j,1);
    int n3 = icone.e(j,2);
    x1.set(&xnod.e(n1-1,0));
    x2.set(&xnod.e(n2-1,0));
    x3.set(&xnod.e(n3-1,0));
    a.set(x2).minus(x1);
    b.set(x3).minus(x1);
    c.cross(a,b);
    double area_elem = 0.5*c.norm_p_all(2.0);
    area_tot += area_elem;
    double area_nod = area_elem/nel;
    for (int k=0; k<nel; k++) {
      int node = icone.e(j,k);
      area.e(node-1) += area_nod;
    }
  }
#if 0
  printf("total area (sum over elems) %f\n",area_tot);
  
  area_tot = 0.0;
  for (int j=0; j<nnod; j++) {
    area_tot += area.e(j);
    // printf("nodo %d, area %f\n",j,area.e(j));
  }
  printf("total area (sum over nodes) %f\n",area_tot);
#endif
}
