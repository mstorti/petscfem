//__INSERT_LICENSE__
// $Id: project4.cpp,v 1.4 2005/03/09 01:52:25 mstorti Exp $

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
#define XNOD1 "square1.nod.tmp"
#define ICONE1 "square1.con.tmp"
#define STATE1 "square1.dat.tmp"
#define XNOD2 "square2.nod.tmp"
#endif

#if 1
#define DATA_DIR "./fluent"
#define XNOD1 DATA_DIR "/fluent.nod"
#define STATE1 DATA_DIR "/fluent.forces"
#define ICONE1 DATA_DIR "/fluent.con"
#define XNOD2 DATA_DIR "/patran.nod"
#endif

#if 0
#define XNOD1  "./mesh1.nod"
#define ICONE1 "./mesh1.con"
#define XNOD2  "./mesh1.nod"
#define STATE1 "./mesh1.nod"
#endif

#if 0
#define XNOD1  "./mesh3.nod"
#define ICONE1 "./mesh3.con"
#define XNOD2  "./mesh3.nod"
#define STATE1 "./mesh3.nod"
#endif

  int ndim = 3;
  int ndimel = 2;
  int nel = 3;
  int ndof = 3;

  dvector<double> xnod1, xnod2, u1, u2,
    area1;
  dvector<int> ico1;

  read_mesh(xnod1,XNOD1, xnod2,XNOD2,
	    u1,STATE1, u2,ico1,ICONE1,
	    ndim,ndimel,nel,ndof);

  nod_vol(xnod1,ico1,area1);
  int nnod = xnod1.size(0);
  for (int j=0; j<nnod; j++)
    for (int k=0; k<ndof; k++)
      u1.e(j,k) /= area1.e(j);

  FemInterp fem_interp;
  fem_interp.init(10,ndof,ndimel,xnod1,ico1);
  u2.clear();
  fem_interp.interp(xnod2,u1,u2);
  u2.print(DATA_DIR "/u2-interp.dat");
}
