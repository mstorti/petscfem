//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.3 2004/11/18 17:44:02 mstorti Exp $

#include <string>

using namespace std;

#include "./femref.h"

int main() { 
  Mesh mesh;
  mesh.read("coord.dat","icone.dat");
}
