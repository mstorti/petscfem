//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.4 2004/11/18 23:34:05 mstorti Exp $

#include <string>

using namespace std;

#include "./femref.h"

int main() { 
  Mesh mesh(2,3);
  mesh.read("coord.dat","icone.dat");
}
