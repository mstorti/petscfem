//__INSERT_LICENSE__
// $Id: project3.cpp,v 1.10 2005/03/09 21:23:23 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <ANN/ANN.h>

#define KNBR 10

// ANNidx nn_idx[KNBR];

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

  FemInterp() :
    kdtree(NULL),
    // nn_idx(NULL),
    // nn_dist(NULL),
    // nn(NULL),
    use_cache(1), tol(1e-6),
    pts(NULL) {}

  void clear() {
    if (kdtree) delete kdtree;
    kdtree = NULL;

#if 0
    if (nn_dist) delete[] nn_dist;
    nn_dist = NULL;
#endif
    nn_dist_v.clear();

#if 0
    if (nn) delete[] nn;
    nn = NULL;
#endif

    // FastMat2::void_cache();

    if (pts) annDeallocPts(pts);
    pts = NULL;

    nn_idx_v.clear();
#if 0
    if (nn_idx) delete[] nn_idx;
    nn_idx = NULL;
#endif
  }

  ~FemInterp() { clear(); }

  void init(int knbr_a, int ndof_a, int ndimel_a,
	    const dvector<double> &xnod_a,
	    const dvector<int> &icone_a);

  void interp(const dvector<double> &xnod2,
	      const dvector<double> &u,
	      dvector<double> &ui);
};

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

  xnod.clone(xnod_a);
  icone.clone(icone_a);

  // Build ANN octree
  FastMat2 xe(1,ndim),xn(1,ndim);
  double inel = 1./nel;
  pts = annAllocPts(nelem,ndim);
  // fixme:= should `pts' be freed after??
  // seems that no
  for (int k=0; k<nelem; k++) {
    xe.set(0.);
    for (int j=0; j<nel; j++) {
      int node = icone.e(k,j);
      xn.set(&xnod.e(node-1,0));
      xe.add(xn);
    }
    xe.scale(inel);
    for (int j=0; j<ndim; j++)
      pts[k][j] = xe.get(j+1);
  }
  kdtree = new ANNkd_tree(pts,nelem,ndim);

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void FemInterp::interp(const dvector<double> &xnod2,
		       const dvector<double> &u,
		       dvector<double> &ui) {
  double tryav=0.0;
  int nnod2 = xnod2.size(0);
  assert(xnod2.size(1) == ndim);
  ui.a_resize(2,nnod2,ndof).defrag();

  double d2min;			// Minimum distance to elements in mesh1
  int k1min;			// Element in mesh1 with minimum distance

  nn_idx_v.resize(knbr);
  ANNidx *nn_idx = &nn_idx_v[0];
  nn_dist_v.resize(knbr);
  ANNdist *nn_dist = &nn_dist_v[0];
  ANNpoint nn = annAllocPt(ndim);

  for (int n2=0; n2<nnod2; n2++) {
    x2.set(&xnod2.e(n2,0));
    if(use_cache) FastMat2::activate_cache(&cache_list);
    for (int j=0; j<ndim; j++)
      nn[j] = xnod2.e(n2,j);
    kdtree->annkSearch(nn,knbr,nn_idx,nn_dist,0.0);
#if 0
    for (int k=0; k<knbr; k++)
      printf("(%d %f) ",nn_idx[k],nn_dist[k]);
    printf("\n");
#endif
    int q;
    for (q=0; q<nelem+knbr; q++) {
      int k = (q<knbr ? nn_idx[q] : q-knbr);
      FastMat2::reset_cache();
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
      assert(ndimel==ndim);

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
      FastMat2::get_cache_position(cp);
      // FastMat2::branch();
      while(true) {
	FastMat2::jump_to(cp);
	// FastMat2::choose(iter);
	iter++;
	invC2.inv(C2);
	L.prod(invC2,b,1,-1,-1);
	// C2.print("C2");
	// L.print("L");
	int neg=0;
	for (int j=1; j<=nel; j++) {
	  double lj = L.get(j);
	  FastMat2::branch();
	  if (lj<-tol) {
	    FastMat2::choose(0);
	    neg=1;
	    int &flag = restricted[j-1];
	    flag = !flag;
	    FastMat2::branch();
	    if (flag) {
	      FastMat2::choose(0);
	      C2.ir(2,j);
	      invCt.ir(2,j);
	      C2.set(invCt).rs();
	      invCt.rs();
	    } else {
	      FastMat2::choose(1);
	      C2.ir(2,j);
	      C.ir(2,j);
	      C2.set(C).rs();
	      C.rs();
	    }
	    FastMat2::leave();
	  }
	  FastMat2::leave();
	}
	if (!neg) break;
	assert(iter<=ndim);
      }
      // FastMat2::leave();
      FastMat2::resync_was_cached();

      // Set area coordinates for restricted
      // nodes to 0.
      for (int j=0; j<nel; j++) {
	FastMat2::branch();
	if (restricted[j]) {
	  FastMat2::choose(0);
	  L.setel(0.,j+1);
	}
	FastMat2::leave();
      }
      for (int j=nel+1; j<=nd1; j++)
	L.setel(0.,j);

      // Set area coordinates for normal
      // vectors to 0.
      C.is(1,1,ndim);
      x2prj.prod(C,L,1,-1,-1);
      C.rs();
      // Form distance vector
      dx2.set(x2).rest(x2prj);
      double d2 = dx2.norm_p_all(2.0);
      if (q==0 || d2<d2min) {
	d2min = d2;
	k1min = k;
	Lmin.set(L);
	x2prjmin.set(x2prj);
      }
      if (d2min<tol) break;
    }
    FastMat2::deactivate_cache();
    // x2.print("x2");
    // x2prjmin.print("x2prjmin");
    // Load state values in `u1_loc'
    for (int j=1; j<=nel; j++) {
      u1_loc.ir(2,j);
      int node = icone.e(k1min,j-1);
      u1_loc.set(&u.e(node-1,0));
    }
    u1_loc.rs();
    // Interpolate
    Lmin.is(1,1,nel);
    u2.prod(u1_loc,Lmin,1,-1,-1);
    Lmin.rs();
    tryav += q+1;
    u2.export_vals(&ui.e(n2,0));
  }
  annDeallocPt(nn);
  // delete[] nn_idx;
  printf("Averg. nbr of tries %f\n",tryav/nnod2);
}

int main() {
#define DATA_DIR "./"
#if 0
#define XNOD1 DATA_DIR "static_p_blade.nod"
#define STATE1 DATA_DIR "static_p_blade.p"
#define ICONE1 DATA_DIR "blade.con"
#define XNOD2 DATA_DIR "patran.nod"
  // #define XNOD2 "./patran.nod"
#elif 0
#define XNOD1 "mesh1.nod"
#define ICONE1 "mesh1.con"
#define STATE1 XNOD1
#define XNOD2 "mesh2.nod"
#elif 0
#define XNOD1 "square1.nod.tmp"
#define ICONE1 "square1.con.tmp"
#define STATE1 "square1.dat.tmp"
#define XNOD2 "square2.nod.tmp"
#else
#define XNOD1 DATA_DIR "fluent.nod"
#define STATE1 DATA_DIR "fluent.forces"
#define ICONE1 DATA_DIR "fluent.con"
#define XNOD2 DATA_DIR "patran.nod"
#endif

  int ndim = 3;
  int ndimel = 3;
  int nel = 3;
  int ndof = 3;

  dvector<double> xnod1,xnod2,u1,u2;
  dvector<int> ico1;

  // Reads mesh1
  xnod1.cat(XNOD1).defrag();
  assert(xnod1.size() % ndim ==0);
  int nnod1 = xnod1.size()/ndim;
  xnod1.reshape(2,nnod1,ndim);
  u1.a_resize(2,nnod1,ndof).read(STATE1);

  ico1.cat(ICONE1).defrag();
  assert(ico1.size() % nel ==0);
  int nelem1 = ico1.size()/nel;
  ico1.reshape(2,nelem1,nel);

  printf("mesh1: %d nodes, %d elems read\n",nnod1,nelem1);

  // Reads mesh2 nodes
  xnod2.cat(XNOD2).defrag();
  assert(xnod2.size() % ndim ==0);
  int nnod2 = xnod2.size()/ndim;
  xnod2.reshape(2,nnod2,ndim);

  printf("mesh2: %d nodes read\n",nnod2);


  FemInterp fem_interp;
  fem_interp.init(KNBR,ndof,ndimel,xnod1,ico1);
  u2.clear();
  fem_interp.interp(xnod2,u1,u2);
  u2.print("field_interpolated.dat");
}
