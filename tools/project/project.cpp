//__INSERT_LICENSE__
// $Id: project.cpp,v 1.1 2005/02/23 23:22:36 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>

int main() {
  dvector<double> x1,x2;
  dvector<int> ico1;
  x1.cat("static_p_blade.nod");
  x1.print();
}
