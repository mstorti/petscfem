//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.4 2005/01/07 20:10:50 mstorti Exp $

#include <string>
#include <list>
#include <multimap.h>
// #include <algorithm>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

struct FaceIterator {
  int elem;
  int face;
};

typedef multimap<int,FaceIterator> face_table_t;

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
  face_table_t face_table;
  GeomObject face, inv_face;
  while (!vis.end()) {  
    GeomObject &go = vis.ref_stack.front().go;
    go.print();
    int nfaces = go.size(GeomObject::OrientedTriT);
    for (int j=0; j<nfaces; j++) {
      go.set(GeomObject::OrientedTriT,j,face);
      face.make_canonical();
      printf("face %d ",j);
      face.print();
      int hash_val = face.hash_val();
      face_table_t::iterator q, 
	q1 = face_table.lower_bound(hash_val),
	q2 = face_table.upper_bound(hash_val);
      if (q1!=q2) {
	printf("possible collition with ");
	for (q = q1; q!=q2; q++)
	  printf("(elem %d, face %d) ", 
		 q->second.elem,q->second.face);
	printf("\n");
      }
    }
    vis.next();
  }
}
