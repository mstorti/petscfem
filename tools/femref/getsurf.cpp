//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.3 2005/01/07 02:39:08 mstorti Exp $

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
  mesh.read("tetra.nod","tetra.con");
#endif
  UniformMesh::visitor vis;
  vis.visit_mode = UniformMesh::BreadthFirst;
  vis.init(mesh);
  while (!vis.end()) {  
    UniformMesh::RefPathNode &w 
      = vis.ref_stack.front();
    w.go.print();
    vis.next();
  }
}
