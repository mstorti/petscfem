//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.1 2005/01/06 18:02:16 mstorti Exp $

#include <string>
#include <list>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

int main() { 

  UniformMesh mesh(OrientedTetraTemplate,3);
#if 0
#define DATA "/u/mstorti/PETSC/GARIBA/DATA"
  mesh.read(DATA "/proy.nod.tmp",DATA "/proy.con.tmp");
#else
  mesh.read("coord.dat","icone.dat");
#endif
  UniformMesh::nat_visitor vis;
  vis.init(mesh);
  while (vis.next()) {  
    UniformMesh::RefPathNode &w 
      = vis.ref_stack.front();
    w.go.print();
  }
}
