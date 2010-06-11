//__INSERT_LICENSE__
// $Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

#include <cstdio>
#include <string>
#include <unistd.h>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <ANN/ANN.h>
#include "./project.h"

extern string print_area_coords;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_mesh(dvector<double> &xnod1,const char *XNOD1,
	       dvector<double> &xnod2,const char *XNOD2,
	       dvector<double> &u1,const char *STATE1,
	       dvector<double> &u2,
	       dvector<int> &ico1,const char *ICONE1,
	       int ndim,int ndimel,int nel,int ndof,
               int use_delaunay) {

  // Reads mesh1
  xnod1.cat(XNOD1).defrag();
  assert(xnod1.size() % ndim ==0);
  int nnod1 = xnod1.size()/ndim;
  xnod1.reshape(2,nnod1,ndim);
  u1.a_resize(2,nnod1,ndof).read(STATE1);

  if (!use_delaunay) {
    ico1.cat(ICONE1).defrag();
    assert(ico1.size() % nel ==0);
    int nelem1 = ico1.size()/nel;
    ico1.reshape(2,nelem1,nel);
    printf("mesh1: %d nodes, %d elems read\n",nnod1,nelem1);
  } else {
    printf("Using delaunay triangulation for interpolation\n");
  }

  // Reads mesh2 nodes
  xnod2.cat(XNOD2).defrag();
  assert(xnod2.size() % ndim ==0);
  int nnod2 = xnod2.size()/ndim;
  xnod2.reshape(2,nnod2,ndim);

  printf("mesh2: %d nodes read\n",nnod2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **argv) {

  MPI_Init(&argc,&argv);
  int ndim = 2;
  int ndimel = 2;
  int nel = ndim+1; // Only for simplices right now
  int ndof = 1;
  int use_delaunay=0;
  print_area_coords = "";

  dvector<double> xnod1, xnod2, u1, u2,
    area1, area2;
  dvector<int> ico1;

  char *xnod1f = strdup("xnod1.tmp");
  char *icone1f = strdup("icone1.tmp");
  char *state1f = strdup("state1.tmp");
  char *xnod2f = strdup("xnod2.tmp");
  char *state2f = strdup("state2.tmp");
  char c;

#define GETOPT_GET(c,fmt,name)			\
    case c:					\
      sscanf(optarg,fmt,&name);			\
      break;
#define SEP "          "
  while ((c = getopt(argc, argv, "hd:l:e:f:x:i:s:y:o:n:a")) != -1) {
    switch (c) {
    case 'h':
      printf(" usage: $ project4.bin -d <NDIM>\n"
	     SEP "-l <NDIMEL> -e <NEL> -f <NDOF>\n"
	     SEP "-x <XNOD1> -i <ICONE1> -s <STATE1>\n"
	     SEP "-y <XNOD2> -o <STATE2> -a\n"
             );
      exit(0);
      GETOPT_GET('d',"%d",ndim);
      GETOPT_GET('l',"%d",ndimel);
      GETOPT_GET('e',"%d",nel);
      GETOPT_GET('f',"%d",ndof);

    case 'x':
      xnod1f = strdup(optarg);
      break;
    case 'i':
      icone1f = strdup(optarg);
      break;
    case 's':
      state1f = strdup(optarg);
      break;
    case 'o':
      state2f = strdup(optarg);
      break;
    case 'y':
      xnod2f = strdup(optarg);
      break;
    case 'n':
      print_area_coords = string(optarg);
      break;
    case 'a':
      use_delaunay = 1;
      break;
    default:
      if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Unknown option character `\\x%x'.\n",
		 optopt);
      abort ();
    }
  }
  PETSCFEM_ASSERT0(!use_delaunay,
                   "Not implemented Delaunay triangulation");  

#if 0
  printf("print_area_coords: %s, size %d\n",
         print_area_coords.c_str(),print_area_coords.size());
#endif

  read_mesh(xnod1,xnod1f, xnod2,xnod2f,
	    u1,state1f,u2,ico1,icone1f,
	    ndim,ndimel,nel,ndof,use_delaunay);

#if 0 // Si los datos vienen `concentrados' por nodos.
  nod_vol(xnod1,ico1,area1);
  int nnod = xnod1.size(0);
  for (int j=0; j<nnod; j++)
    for (int k=0; k<ndof; k++)
      u1.e(j,k) /= area1.e(j);
#endif

  FemInterp fem_interp;
  if (!use_delaunay) {
    fem_interp.init(10,ndof,ndimel,xnod1,ico1);
  } else {
    assert(!strcmp(icone1f,"icone1.tmp"));
    assert(nel==ndim+1);
    fem_interp.init(10,ndof,ndimel,xnod1);
  }
  u2.clear();
  fem_interp.interp(xnod2,u1,u2);
  u2.print(state2f);

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
  MPI_Finalize();
  return 0;
}
