//__INSERT_LICENSE__
// $Id mstorti-v6-2-1-g23d6622 Wed Jun 20 11:58:02 2007 -0300$

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <ANN/ANN.h>
#include "./project.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_mesh(dvector<double> &xnod1,const char *XNOD1,
	       dvector<double> &xnod2,const char *XNOD2,
	       dvector<double> &u1,const char *STATE1,
	       dvector<double> &u2,
	       dvector<int> &ico1,const char *ICONE1,
	       int ndim,int ndimel,int nel,int ndof) {

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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main() {

#if 0
#define DATA_DIR1 "/home/mstorti/PETSC/petscfem-cases/sqcav-ther-Ra1.6e9-N100"
#define DATA_DIR "/home/mstorti/PETSC/petscfem-cases/sqcav-ther"
#define XNOD1  DATA_DIR1 "/sqcav-ther.nod.tmp"
#define ICONE1 DATA_DIR1 "/sqcav-ther.con-tri.tmp"
#define STATE1 DATA_DIR "/sqcav-ther.state-7418"
#define XNOD2  DATA_DIR "/sqcav-ther.nod.tmp"
#endif

#if 1
#define XNOD1  "./xnod1.tmp"
#define ICONE1 "./icone1.tmp"
#define STATE1 "./u.tmp"
#define XNOD2  "./xnod2.tmp"
#endif

  int ndim = 2;
  int ndimel = 2;
  int nel = ndim+1; // Only for simplices right now
  int ndof = 1;

  dvector<double> xnod1, xnod2, u1, u2,
    area1, area2;
  dvector<int> ico1;

  read_mesh(xnod1,XNOD1, xnod2,XNOD2,
	    u1,STATE1, u2,ico1,ICONE1,
	    ndim,ndimel,nel,ndof);

#if 0 // Si los datos vienen `concentrados' por nodos.
  nod_vol(xnod1,ico1,area1);
  int nnod = xnod1.size(0);
  for (int j=0; j<nnod; j++)
    for (int k=0; k<ndof; k++)
      u1.e(j,k) /= area1.e(j);
#endif

  FemInterp fem_interp;
  fem_interp.init(10,ndof,ndimel,xnod1,ico1);
  u2.clear();
  fem_interp.interp(xnod2,u1,u2);
  u2.print("./u2-interp.dat");

#if 0
  dvector<int> ico2;
  ico2.cat(ICONE2).defrag();
  assert(ico2.size() % nel ==0);
  int nelem2 = ico2.size()/nel;
  ico2.reshape(2,nelem2,nel);
  printf("mesh 2, %d elements read\n",nelem2);

  nod_vol(xnod2,ico2,area2);
  int nnod2 = xnod2.size(0);
  for (int j=0; j<nnod2; j++)
    for (int k=0; k<ndof; k++)
      u2.e(j,k) *= area2.e(j);
  u2.print("u2n.dat");
#endif
}
