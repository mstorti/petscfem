/*__INSERT_LICENSE__*/
// $Id: test10.cpp,v 1.3 2003/08/10 13:54:14 mstorti Exp $

#include <cstdio>
#include <src/dvector.h>
#include <src/dvector2.h>

int main(int argc, char **argv) {
  dvector<double> v,w,x;
  v.cat("data.txt.tmp").print("data.out.tmp");
  w.mono(1000,0.).set(2.);
  x.mono(2000,0.).set(3.);
  printf("ok\n");
}
