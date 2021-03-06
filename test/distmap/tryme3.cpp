//__INSERT_LICENSE__
//$Id: tryme3.cpp,v 1.2 2003/01/08 15:49:04 mstorti Exp $
#include <iostream>

class A {
  int *a;
public:
  A() {a = new int[1000];};
  virtual ~A() {delete a; cout << "~A() \n";};
};

class B : public A {
  int *b;
public:
  B() {b = new int[1000];};
  ~B() {delete b; cout << "~B() \n";};
};

main() {
  A *a;
  for (int j=0; j<10; j++) {
    a = new B;
    // The question is if `delete a' calls the destructor
    // of B or not 
    delete a;
  }
}
