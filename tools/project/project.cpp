//__INSERT_LICENSE__
// $Id: project.cpp,v 1.19 2005/02/25 20:04:33 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>

static double drand() { 
  return double(rand())/double(RAND_MAX); 
}

int main() {
  int ndim = 3;
  int nel = 3;
  int ndof = 1;

  dvector<double> xnod1,xnod2,u1;
  dvector<int> ico1;

#define DATA_DIR "/u/mstorti/PETSC/tatuus/data/"
#if 1
#define XNOD1 DATA_DIR "static_p_blade.nod"
#define STATE1 DATA_DIR "static_p_blade.p"
#define ICONE1 DATA_DIR "blade.con"
#define XNOD2 DATA_DIR "patran.nod"
  // #define XNOD2 "./patran.nod"
#else
#define XNOD1 "mesh1.nod"
#define ICONE1 "mesh1.con"
#define XNOD2 "mesh2.nod"
#endif

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

  // coordinates matrix. 
  int nd1 = ndim+1;
  FastMat2 C(2,nd1,nd1),C2(2,nd1,nd1),
    invC(2,nd1,nd1), invC2(2,nd1,nd1),
    invCt(2,nd1,nd1),
    x2(1,ndim),dx2(1,ndim),
    x2prj(1,ndim),x2prjmin(1,ndim), x12(1,ndim),
    x13(1,ndim),x1(1,ndim),
    nor(1,ndim),L(1,nd1),
    b(1,ndim+1),u1_loc(2,ndof,nel),
    u2(1,ndof);
  vector<int> restricted(ndim);
  double tol = 1e-6;
  double d2min;			// Minimum distance to elements in mesh1
  int k1min;			// Element in mesh1 with minimum distance
  FastMat2 Lmin(1,nd1);
  FastMatCachePosition cp,cp1;
  FastMatCacheList cache_list;
  int use_cache;
  FILE *fid = fopen("pinterp.dat","w");
#if 0
  printf("use cache > ");
  scanf("%d",&use_cache);
  printf("use_cache %d\n",use_cache);
#else
  use_cache = 1;
#endif
  for (int n2=0; n2<nnod2; n2++) {
    x2.set(&xnod2.e(n2,0));
    if(use_cache) FastMat2::activate_cache(&cache_list);
    for (int k=0; k<nelem1; k++) {
      FastMat2::reset_cache();
      C.is(1,1,ndim);
      for (int j=0; j<nel; j++) {
	int node = ico1.e(k,j);
	C.ir(2,j+1).set(&xnod1.e(node-1,0));
	if (j==0) {
	  x12.set(C);
	  x13.set(C);
	} else if (j==1) x12.rest(C);
	else if (j==2) x13.rest(C);
      }
      x12.scale(-1.0);
      x13.scale(-1.0);
      nor.cross(x12,x13);
      C.rs();
      C.ir(1,nd1).is(1,1,ndim).set(1.0).rs();
      C.ir(2,nd1).is(1,1,ndim).set(nor).rs();
      C.setel(0.0,nd1,nd1);
      b.is(1,1,ndim).set(x2).rs();
      b.setel(1.0,nd1);
      for (int j=0; j<ndim; j++)
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
	for (int j=1; j<=ndim; j++) {
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
      // printf("converged on iters %d\n",iter);
      for (int j=0; j<ndim; j++) {
	FastMat2::branch();
	if (restricted[j]) {
	  FastMat2::choose(0);
	  L.setel(0.,j+1);
	}
	FastMat2::leave();
      }
      L.setel(0.,nd1);
      C.is(1,1,ndim);
      x2prj.prod(C,L,1,-1,-1);
      C.rs();
      // Form distance vector
      dx2.set(x2).rest(x2prj);
      double d2 = dx2.norm_p_all(2.0);
      if (k==0 || d2<d2min) {
	d2min = d2;
	k1min = k;
	Lmin.set(L);
	x2prjmin.set(x2prj);
      }
    }
    FastMat2::deactivate_cache();
    // x2.print("x2");
    // x2prjmin.print("x2prjmin");
    // Load state values in `u1_loc'
    for (int j=1; j<=nel; j++) {
      u1_loc.ir(2,j);
      int node = ico1.e(k1min,j-1);
      u1_loc.set(u1.e(node-1,0));
    }
    u1_loc.rs();
    // Interpolate
    Lmin.is(1,1,ndim);
    u2.prod(u1_loc,Lmin,1,-1,-1);
    Lmin.rs();
    printf("node2 %d, dist min %f\n",n2,d2min);
    for (int j=1; j<=ndof; j++)
      fprintf(fid,"%f ",u2.get(j));
    fprintf(fid,"\n");
    u2.print("u2");
  }
  FastMat2::void_cache();
  fclose(fid);
}
