/*__INSERT_LICENSE__*/
// $Id: sp.cpp,v 1.1 2001/09/20 20:13:28 mstorti Exp $

#include <sparse.h>

using namespace Sparse;

int main() {
  Vec v(5);
  v.set(2,1.).set(3,2.).print();
  v.set(10,10.).set(11,11.).print();
  v.grow(0).set(5,5.).print();
  v.set(20,20.).print();
  
}
