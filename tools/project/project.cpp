//__INSERT_LICENSE__
// $Id: project.cpp,v 1.7 2005/02/25 00:28:13 mstorti Exp $

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

  // Constructs `Sinv_v' a vector of doubles storing
  // the coefficients of the inverse of element
  // coordinates matrix. 
  int nd1 = ndim+1;
  FastMat2 C(2,nd1,nd1),invC(2,nd1,nd1);
  for (int n2=0; n2<
  for (int k=0; k<nelem1; k++) {
    // Nodes of element (base 0)
    for (int j=0; j<nel; j++) {
      int node = ico1.e(k.
      n1 = ico1.e(k,0),
      n2 = ico1.e(k,1),
      n3 = ico1.e(k,2);
    // Position of first node
    xn1.set(&xnod1.e(n1-1,0)) ;
    xn2.set(&xnod1.e(n2-1,0)) ;
    xn3.set(&xnod1.e(n3-1,0)) ;
    x12.set(xn2).rest(xn1);
    S.ir(2,2).set(x12);
    x13.set(xn3).rest(xn1);
    S.ir(2,3).set(x13);
    S.ir(2,1).cross(x12,x13).rs();
    Sinv.inv(S);
    Sinv.export_vals(&Sinv_v.e(k,0,0));
  }


  FastMat2 x2(1,ndim),L(1,ndim),
    dx(1,ndim),xcloser(1,ndim);
  // Loop over mesh2 nodes
  for (int k2=0; k2<nnod2; k2++) {
    // A node in the second mesh (test node)
    x2.set(&xnod2.e(k2,0));
    // Loop over elements in the first mesh
    for (int e1=0; e1<nelem1; e1++) {
      // First note of the element
      int n1 = ico1.e(e1,0);
      // Coordinate of the first node
      xn1.set(&xnod1.e(n1-1,0));
      // Coordinate of the test node w.r.t.
      // first node of element
      dx.set(x2).rest(xn1);
      // inverse basis for the element
      Sinv.set(&Sinv_v.e(e1,0,0));
      // Area coordinates
      L.prod(Sinv,dx,1,-1,-1);
      double *lv = L.storage_begin();
      // Redefine third area coordinate
      lv[0] = 1.0-lv[1]-lv[2];
      printf("x %f %f %f, L %f %f %f\n",
	     x2.get(1),x2.get(2),x2.get(3),
	     L.get(1),L.get(2),L.get(3));
      // Number of negative area coordinates
      int neg=0, // nbr of negative area coordinates
	pos=-1; //some area coordinate that is positive
      for (int j=0; j<ndim; j++) {
	neg += lv[j]<0;
	if (lv[j]>=0) pos = j;
      }
      assert(neg<ndim);
      assert(pos>=0);
      if (neg==0) {
	xcloser.prod(S,
      } else if (neg==1) {
	int indx;
	if (lv[0]<0) indx=1;
	else if (lv[1]<0) indx=2;
	else if (lv[2]<0) indx=3;
	
      }
    }
  }
}
