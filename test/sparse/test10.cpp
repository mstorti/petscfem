/*__INSERT_LICENSE__*/
// $Id: test10.cpp,v 1.1 2003/08/09 13:42:25 mstorti Exp $

#include <src/dvector.h>
#include <src/dvector2.h>

int main(int argc, char **argv) {
  dvector<double> v;
  v.cat("data.txt.tmp").print("data.out.tmp");
}
