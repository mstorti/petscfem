/*__INSERT_LICENSE__*/
// $Id: test11.cpp,v 1.1 2006/04/08 20:36:13 mstorti Exp $

#include <cstdio>

#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/dvecpar.h>

int main(int argc, char **argv) {
  PetscInitialize(&argc,&argv,NULL,NULL);
  dvector<double> v;
  dvector_read_parallel("file.dat",v);
  double sum = 0.0;
  for (int k=0; k<v.size(); k++) 
    sum += v.ref(k);
  printf("sum %f\n",sum);
  PetscFinalize();
}
