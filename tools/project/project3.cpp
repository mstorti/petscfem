//__INSERT_LICENSE__
// $Id: project3.cpp,v 1.9 2005/03/02 23:26:31 mstorti Exp $

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

#define XNOD1 "square1.nod.tmp"
#define ICONE1 "square1.con.tmp"
#define STATE1 "square1.dat.tmp"
#define XNOD2 "square2.nod.tmp"

  int ndim = 2;
  int ndimel = 2;
  int nel = 3;
  int ndof = 2;

  dvector<double> xnod1, xnod2, u1, u2;
  dvector<int> ico1;

  read_mesh(xnod1,"square1.nod.tmp",
	    xnod2,"square2.nod.tmp",
	    u1,"square1.dat.tmp",
	    u2,ico1,"square1.con.tmp",
	    ndim,ndimel,nel,ndof);
  
  dvector<double> xnod1b, xnod2b, u1b, u2b;
  dvector<int> ico1b;

  read_mesh(xnod1b,"squareb1.nod.tmp",
	    xnod2b,"squareb2.nod.tmp",
	    u1b,"squareb1.dat.tmp",
	    u2b,ico1b,"squareb1.con.tmp",
	    ndim,ndimel,nel,ndof);

#define KNBR 10
  while (1) {
    FemInterp fem_interp;
    fem_interp.init(KNBR,2,2,xnod1,ico1);
    u2.clear();
    fem_interp.interp(xnod2,u1,u2);

    FemInterp fem_interpb;
    fem_interpb.init(KNBR,2,2,xnod1b,ico1b);
    u2b.clear();
    fem_interpb.interp(xnod2b,u1b,u2b);
  }
}
