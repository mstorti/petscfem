/*__INSERT_LICENSE__*/
//$Id: testfm2c.cpp,v 1.3.8.1 2001/12/21 00:13:32 mstorti Exp $
 
//< if (0) { //>//
#include <stdio.h>
#include <time.h>

#include <sles.h>

#include <src/fastmat2.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>

#define XNOD(j,k) (xnod[2*(j)+(k)])
#define ICONE(j,k) (icone[3*(j)+(k)])

// Test and example of use for the FastMat2 matrix pacakage.
// Compute the minimu size of a FEM mesh of triangular elements.

int main() {

  int n=300; // Number of elements por side of the squares, so that
	    // we have 2*n^2 elements
  double noise=0.3; // Amount of noise added to the mesh

  int nnod=(n+1)*(n+1);
  int nelem=2*n*n;
  double hav = 1./double(n); //typical size of the element

  double *xnod = new double[2*nnod];
  int *icone = new int[3*nelem];

  // define position of nodes
  int node=0;
  for (int k=0; k<n+1; k++) {
    double y=double(k)/double(n);
    if (k>0 && k<n) y += noise*(2*drand()-1)*hav;
    for (int j=0; j<n+1; j++) {
      double x=double(j)/double(n);
      if (j>0 && j<n) x += noise*(2*drand()-1)*hav;
      XNOD(node,0)=x;
      XNOD(node,1)=y;
      node++;
    }
  }

  int iele=0;
  double minarea,maxarea,total_area=0;
  for (int k=0; k<n; k++) {
    for (int j=0; j<n; j++) {
      int n1,n2,n3,n4;
      n1=k*(n+1)+j+1;
      n2=n1+1;
      n3=n2+n+1;
      n4=n1+n+1;

      ICONE(iele,0)=n1;
      ICONE(iele,1)=n2;
      ICONE(iele,2)=n3;
      iele++;

      ICONE(iele,0)=n3;
      ICONE(iele,1)=n4;
      ICONE(iele,2)=n1;
      iele++;
    }
  }
//< } //>//
  Chrono chrono;
  FastMat2 x(2,3,2),a(1,2),b(1,2),J(2,2,2);
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  // Compute area of elements
  chrono.start();
  for (int ie=0; ie<nelem; ie++) {
    FastMat2::reset_cache();
    for (int k=1; k<=3; k++) {
      int node = ICONE(ie,k-1);
      x.ir(1,k).set(&XNOD(node-1,0)).rs();
    }
    x.rs();
    a.set(x.ir(1,2));
    a.rest(x.ir(1,1));

    b.set(x.ir(1,3));
    b.rest(x.ir(1,1));

    J.ir(1,1).set(a);
    J.ir(1,2).set(b).rs();
    
    double area = J.det()/2.;
    total_area += area;
    if (ie==0) {
      minarea = area;
      maxarea = area;
    }

    if (area>maxarea) maxarea=area;
    if (area<minarea) minarea=area;
  }
  printf("total_area %g, min area %g,max area %g, ratio: %g\n",
	 total_area,minarea,maxarea,maxarea/minarea);
  printf("Total area OK? : %s\n",
	 (fabs(total_area-1)<1e-8 ? "YES" : "NOT"));
  double cpu = chrono.elapsed();
  FastMat2::print_count_statistics();
  printf("CPU: %g, number of elements: %d\n"
	 "rate: %g [sec/Me], %g Mflops\n",
	 cpu,nelem,cpu*1e6/nelem,
	 nelem*FastMat2::operation_count()/cpu/1e6);
  FastMat2::deactivate_cache();
  FastMat2::void_cache();
//< if (0) { //>//
}
//< } //>//
