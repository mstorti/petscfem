//__INSERT_LICENSE__
// $Id: project.cpp,v 1.10 2005/02/25 01:16:13 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>

int main() {
  int ndim = 3;
  int nel = 3;

  dvector<double> xnod1,xnod2;
  dvector<int> ico1;

#define DATA_DIR "/home/mstorti/PETSC/tatuus/data/"
#if 0
#define XNOD1 DATA_DIR "static_p_blade.nod"
#define ICONE1 DATA_DIR "blade.con"
#define XNOD2 DATA_DIR "static_p_blade.nod"
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
    x2(1,ndim),x12(1,ndim),
    x13(1,ndim),x1(1,ndim),
    nor(1,ndim),L(1,ndim+1),
    b(1,ndim+1);
  vector<int> restricted(ndim);
  for (int n2=0; n2<nnod2; n2++) {
    x2.set(&xnod2.e(n2,0));
    for (int k=0; k<nelem1; k++) {
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
      b.print("b");
      C.print("C");
      invCt.ctr(invC,2,1).scale(-1.0);
      invCt.ir(1,nd1).set(0.).rs();
      invCt.print("invCt");
      int iter=0;
      while(true) {
	iter++;
	invC2.inv(C2);
	L.prod(invC2,b,1,-1,-1);
	C2.print("C2");
	L.print("L");
	int neg=0;
	for (int j=1; j<=ndim; j++) {
	  if (L.get(j)<0) {
	    neg=1;
	    restricted[j-1]=1;
	    C2.ir(2,j);
	    invCt.ir(2,j);
	    C2.set(invCt);
	  }
	}
	if (!neg || iter>ndim) break;
	C2.rs();
	invCt.rs();
      }
      printf("converged on iters %d\n",iter);
    }
  }
}
