/*__INSERT_LICENSE__*/
// $Id: test10.cpp,v 1.2 2003/08/10 01:32:06 mstorti Exp $

#include <cstdio>
#include <src/dvector.h>
#include <src/dvector2.h>

int main(int argc, char **argv) {
  dvector<double> v,w;
  v.cat("data.txt.tmp").print("data.out.tmp");
  w.mono(1000,0.);
  printf("ok\n");
}
