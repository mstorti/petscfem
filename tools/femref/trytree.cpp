//__INSERT_LICENSE__
// $Id: trytree.cpp,v 1.1 2004/12/26 00:41:00 mstorti Exp $

#include <cstdlib>
#include <cassert>
#include <iostream>
#include "./tree.h"

using namespace std;

int main() {
  tree<int> A;
  typedef tree<int>::iterator node_t;
  node_t n = A.begin();
  n = A.insert(n,10);
  A.lisp_print();
}
