/*__INSERT_LICENSE__*/
//$Id: testfm2d.cpp,v 1.3 2001/05/30 03:58:53 mstorti Exp $

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#include <sles.h>

#include "../src/fastmat2.h"
#include "../src/fastmat.h"
#include "../src/fem.h"
#include "../src/utils.h"
#include "../src/util2.h"
#include <newmatio.h>

double gettod() 
{
 struct timeval tv;
 gettimeofday(&tv,0);
 return tv.tv_sec + 1e-6 * tv.tv_usec;
}

void rand_init(FastMat2 & A) {
  Indx dims;
  A.get_dims(dims);
  int n = mem_size(dims);
  double * r = new double[n];
  for (int j=0; j<n; j++) {
    r[j] = 2.*drand()-1.;
  }
  A.set(r);
  delete[] r;
}

inline double * location(double *a, int i, int j, int m) {
  return a + m * i + j;
}

void myprod(double *c, double *a, double *b,int m, int p, int n) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      double *pc = location(c,i,j,n);
      double *pa = location(a,i,0,p);
      double *pb = location(b,0,j,n);
      double sum = 0.;
      double *paend = pa + n;
      for (int k=0; k<p; k++) {
	//      while (pa<paend) {
	sum += *pa * *pb;
	pa++;
	pb += n;
      }
      *pc = sum;
    }
  }
}

void myprodt(double *c, double *a, double *b,int m, int p, int n) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      double *pc = location(c,i,j,n);
      double *pa = location(a,i,0,p);
      double *pb = location(b,j,0,n);
      double sum = 0.;
      double *paend = pa + n;
      for (int k=0; k<p; k++) {
	// while (pa<paend) {
	sum += *pa++ * *pb++;
      }
      *pc = sum;
    }
  }
}

void something(FastMat2 & A) {
  A.setel(23.,1,1);
}

int main (int argc, char **argv) {
  Chrono chrono;
  FastMatCacheList cache_list;
  
  int n=30,N=100,M,c,Mset=0,fm2=1,t=5;
  M=N;
  while ((c = getopt (argc, argv, "n:N:M:Ft:")) != -1) {
    switch (c) {
    case 'F':
      fm2=0;
      break;
    case 'n':
      sscanf(optarg,"%d",&n);
      break;
    case 'N':
      sscanf(optarg,"%d",&N);
      break;
    case 'M':
      sscanf(optarg,"%d",&M);
      Mset=1;
      break;
    case 't':
      sscanf(optarg,"%d",&t);
      break;
    default:
      abort ();
    }
  }

  if (!Mset) M=N;

  // int n=30, N=10000;
  FastMat2 A(2,n,n),C,D(2,n,n);
  FastMat AA(n,n),CC(n,n);
  int NN = M*n*n; 
  double opstyp=double(t)*1e8;
  // double opstyp=30.*30.*30.*100.*100.;
  int nloop = int(opstyp/(n*n*n*N));
  printf("n: %d, N: %d, M: %d, loop: %d\n",n,N,M,nloop);
  
  double *a=new double[NN];
  for (int j=0; j<NN; j++) 
    a[j] = 2*drand()-1;

  printf("fm2: %d\n",fm2);
  int cpy_it_max = 3;
  double *cpu_v = new double[cpy_it_max];
  for (int cpy_it=1; cpy_it<=cpy_it_max; cpy_it++) {
    chrono.start();
    double start=gettod();
    for (int jl=0; jl<nloop; jl++) {
      if (fm2) {
	FastMat2::activate_cache(&cache_list);
	for (int j=0; j<N; j++) {
	  FastMat2::reset_cache();
	  for (int kkk=0; kkk<cpy_it; kkk++) {
	    A.set(&a[n*n*(j % M)]);
	    //A.set(a);
	  }
	  C.prod(A,A,1,-1,-1,2);
//  	  A.print("A: ");
//  	  C.print("C: ");
//  	  exit(0);
	}
	FastMat2::void_cache();
	FastMat2::deactivate_cache();
      } else {
	assert(0); // Not valid right now !!
	for (int j=0; j<N; j++) {
	  // AA.set(a);
	  A.set(&a[n*n*(j % M)]);
	}
	myprodt(CC.store,AA.store,AA.store,n,n,n);
      }
    }

    double delta=gettod()-start;
    double cpu = chrono.elapsed();
    printf("cpy times: %d, cpu: %f\n",cpy_it,cpu);
    cpu_v[cpy_it-1] = cpu;
    if (fabs(delta-cpu)/cpu>.1) 
      printf("warning: elapsed and cpu differ too much!!\n"
	     "elapsed: %f, cpu: %f [Mflops]\n",delta,cpu);
  }

  double cpu = 2*cpu_v[0]-cpu_v[1];
  for (int cpy_it=3; cpy_it<=cpy_it_max; cpy_it++) {
    double cpuu = cpu_v[0]-(cpu_v[cpy_it-1]-cpu_v[0])/double(cpy_it-1);
    double dev = fabs(cpuu-cpu)/cpu;
    if (dev>0.1) 
      printf("warning: deviation from linear exceeds tol., dev= %f\n",dev);
  }
  double ops = (double) (2.*nloop*N*n*n*n);
  printf("elapsed: %f [sec], ops: %f [flop], rate: %f [Mflops]\n",cpu,ops,ops/cpu/1e6);
  
}
