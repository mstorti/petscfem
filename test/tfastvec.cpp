//__INSERT_LICENSE__
//$Id: tfastvec.cpp,v 1.2 2001/04/01 01:35:07 mstorti Exp $
 
#include <iostream>
#include <stdio.h>
#include "../src/fastlib.h"

int main() {
  // dimension a vector of 10 components
  int n=10;
  FastVector<int> a(n,0);
  for (int j=0; j<n; j++) {
    a[j] = j;
  }

  // resizes to 40
  int nn=40;
  a.resize(nn);
  for (int j=0; j<nn; j++) {
    a[j] = j;
  }
  a.print("resized");

  // resizes back to 20
  nn=20;
  a.resize(nn);
  for (int j=0; j<nn; j++) {
    a[j] = j;
  }
  a.print("resized");

  FastVector<double> c(n,3.3);
  c.print("with doubles");

}
